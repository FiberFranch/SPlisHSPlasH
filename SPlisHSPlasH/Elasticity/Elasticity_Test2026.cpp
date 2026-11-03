// #define ELASTICITY_DEBUG_SINGLE_THREAD
#include "Elasticity_Test2026.h"
#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "Utilities/Timing.h"

using namespace SPH;
using namespace GenParam;

std::string Elasticity_Test2026::METHOD_NAME = "Elasticity Test 2026";
int Elasticity_Test2026::YOUNGS_MODULUS = -1;
int Elasticity_Test2026::POISSON_RATIO = -1;
int Elasticity_Test2026::FIXED_BOX_MIN = -1;
int Elasticity_Test2026::FIXED_BOX_MAX = -1;

Elasticity_Test2026::Elasticity_Test2026(FluidModel *model) : NonPressureForceBase(model)
{
    const unsigned int numParticles = model->numActiveParticles();
    m_stress.resize(numParticles);

    m_youngsModulus = static_cast<Real>(100000.0);
    m_poissonRatio = static_cast<Real>(0.3);
    m_fixedBoxMin.setZero();
    m_fixedBoxMax.setZero();

    model->addField({ "deviatoric stress", METHOD_NAME, FieldType::Vector6, [&](const unsigned int i) -> Real* { return &m_stress[i][0]; }, true });
}

Elasticity_Test2026::~Elasticity_Test2026(void)
{
    m_model->removeFieldByName("deviatoric stress");
}

void Elasticity_Test2026::deferredInit()
{
    initValues();
}

void Elasticity_Test2026::initParameters()
{
    NonPressureForceBase::initParameters();

    YOUNGS_MODULUS = createNumericParameter("youngsModulus", "Young`s modulus", &m_youngsModulus);
    setGroup(YOUNGS_MODULUS, "Fluid Model|Elasticity");
    setDescription(YOUNGS_MODULUS, "Stiffness of the elastic material");
    RealParameter* rparam = static_cast<RealParameter*>(getParameter(YOUNGS_MODULUS));
    rparam->setMinValue(0.0);

    POISSON_RATIO = createNumericParameter("poissonsRatio", "Poisson`s ratio", &m_poissonRatio);
    setGroup(POISSON_RATIO, "Fluid Model|Elasticity");
    setDescription(POISSON_RATIO, "Ratio of transversal expansion and axial compression");
    rparam = static_cast<RealParameter*>(getParameter(POISSON_RATIO));
    rparam->setMinValue(static_cast<Real>(-1.0 + 1e-4));
    rparam->setMaxValue(static_cast<Real>(0.5 - 1e-4));

    ParameterBase::GetVecFunc<Real> getFctFMin = [&]()-> Real* { return m_fixedBoxMin.data(); };
    ParameterBase::SetVecFunc<Real> setFctFMin = [&](Real* val)
    {
        m_fixedBoxMin = Vector3r(val[0], val[1], val[2]);
        determineFixedParticles();
    };
    FIXED_BOX_MIN = createVectorParameter("fixedBoxMin", "Fixed box min", 3u, getFctFMin, setFctFMin);
    setGroup(FIXED_BOX_MIN, "Fluid Model|Elasticity");
    setDescription(FIXED_BOX_MIN, "Minimum point of box of which contains the fixed particles.");
    getParameter(FIXED_BOX_MIN)->setReadOnly(true);

    ParameterBase::GetVecFunc<Real> getFctFMax = [&]()-> Real* { return m_fixedBoxMax.data(); };
    ParameterBase::SetVecFunc<Real> setFctFMax = [&](Real* val)
    {
        m_fixedBoxMax = Vector3r(val[0], val[1], val[2]);
        determineFixedParticles();
    };
    FIXED_BOX_MAX = createVectorParameter("fixedBoxMax", "Fixed box max", 3u, getFctFMax, setFctFMax);
    setGroup(FIXED_BOX_MAX, "Fluid Model|Elasticity");
    setDescription(FIXED_BOX_MAX, "Maximum point of box of which contains the fixed particles.");
    getParameter(FIXED_BOX_MAX)->setReadOnly(true);

}

/** Mark all particles in the bounding box as fixed.
*/
void Elasticity_Test2026::determineFixedParticles()
{
    const unsigned int numParticles = m_model->numActiveParticles();
    int numFixed = 0;

    if (!m_fixedBoxMin.isZero() || !m_fixedBoxMax.isZero())
    {
        for (int i = 0; i < (int)numParticles; i++)
        {
            const Vector3r& x = m_model->getPosition0(i);
            if ((x[0] > m_fixedBoxMin[0]) && (x[1] > m_fixedBoxMin[1]) && (x[2] > m_fixedBoxMin[2]) &&
                (x[0] < m_fixedBoxMax[0]) && (x[1] < m_fixedBoxMax[1]) && (x[2] < m_fixedBoxMax[2]))
            {
                m_model->setParticleState(i, ParticleState::Fixed);
                numFixed++;
            }
        }
    }
    LOG_INFO << "Elasticity_Test2026: Marked " << numFixed << " particles as fixed.";
}

void Elasticity_Test2026::initValues()
{
    determineFixedParticles();
    findObjects();
}

void Elasticity_Test2026::findObjects()
{
    const unsigned int numParticles = m_model->numActiveParticles();
    m_objects.clear();
    std::map<unsigned int, unsigned int> obj2group;

    for (unsigned int i = 0; i < numParticles; i++)
    {
        const unsigned int objId = m_model->getObjectId(i);

        if (obj2group.find(objId) == obj2group.end())
        {
            m_objects.push_back(new ElasticObject());
            const unsigned int groupIndex = (unsigned int)(m_objects.size() - 1);
            obj2group[objId] = groupIndex;
            m_objects[groupIndex]->m_particleIndices.push_back(i);
        }
        else
        {
            const unsigned int groupIndex = obj2group[objId];
            m_objects[groupIndex]->m_particleIndices.push_back(i);
        }
    }

    for (size_t groupIndex = 0; groupIndex < m_objects.size(); groupIndex++)
    {
        struct Comparator {
            Comparator(Elasticity_Test2026* _this) : m_this(_this) {};
            bool operator()(unsigned int a, unsigned int b)
            {
                if ((m_this->m_model->getParticleState(a) != ParticleState::Active) && (m_this->m_model->getParticleState(b) == ParticleState::Active))
                    return false;
                else if ((m_this->m_model->getParticleState(a) == ParticleState::Active) && (m_this->m_model->getParticleState(b) != ParticleState::Active))
                    return true;
                else
                    return a < b;
            }
            Elasticity_Test2026* m_this;
        };

        std::sort(m_objects[groupIndex]->m_particleIndices.begin(), m_objects[groupIndex]->m_particleIndices.end(), Comparator(this));

        m_objects[groupIndex]->m_nFixed = 0;
        for (size_t i = 0; i < m_objects[groupIndex]->m_particleIndices.size(); i++)
        {
            if (m_model->getParticleState(m_objects[groupIndex]->m_particleIndices[i]) != ParticleState::Active)
                m_objects[groupIndex]->m_nFixed++;
        }
        LOG_INFO << "Object " << groupIndex << " - fixed particles: " << m_objects[groupIndex]->m_nFixed;
    }
}

void Elasticity_Test2026::step()
{
    computeStress();
    computeForces();
}

void Elasticity_Test2026::computeStress()
{
    const unsigned int numParticles = m_model->numActiveParticles();
    if (numParticles == 0) return;

    const Real dt = TimeManager::getCurrent()->getTimeStepSize();

    Simulation *sim = Simulation::getCurrent();
    const unsigned int fluidModelIndex = m_model->getPointSetIndex();

    Real mu = m_youngsModulus / (static_cast<Real>(2.0) * (static_cast<Real>(1.0) + m_poissonRatio));

#ifndef ELASTICITY_DEBUG_SINGLE_THREAD
    #pragma omp parallel default(shared)
    {
        #pragma omp for schedule(static)
#endif
        for (int i = 0; i < (int)numParticles; i++)
        {
            if (m_model->getParticleState(i) == ParticleState::Active)
            {
                const Vector3r &xi = m_model->getPosition(i);
                const Vector3r &vi = m_model->getVelocity(i);
                Matrix3r gradV = Matrix3r::Zero();

                const unsigned int numNeighbors = sim->numberOfNeighbors(fluidModelIndex, fluidModelIndex, i);

                for (unsigned int j = 0; j < numNeighbors; j++)
                {
                    const unsigned int neighborIndex = sim->getNeighbor(fluidModelIndex, fluidModelIndex, i, j);
                    const Vector3r &xj = m_model->getPosition(neighborIndex);
                    const Vector3r &vj = m_model->getVelocity(neighborIndex);
                    const Real Vj = m_model->getMass(neighborIndex) / m_model->getDensity(neighborIndex);

                    gradV += Vj * (vj - vi) * sim->gradW(xi - xj).transpose();
                }

                Matrix3r strainRate = static_cast<Real>(0.5) * (gradV + gradV.transpose());
                Matrix3r spinRate = static_cast<Real>(0.5) * (gradV - gradV.transpose());

                Matrix3r stressMat;
                stressMat(0, 0) = m_stress[i][0];
                stressMat(1, 1) = m_stress[i][1];
                stressMat(2, 2) = m_stress[i][2];
                stressMat(0, 1) = stressMat(1, 0) = m_stress[i][3];
                stressMat(0, 2) = stressMat(2, 0) = m_stress[i][4];
                stressMat(1, 2) = stressMat(2, 1) = m_stress[i][5];

                Matrix3r stressRate = static_cast<Real>(2.0) * mu * strainRate;
                Matrix3r jaumannRate = stressRate + stressMat * spinRate - spinRate * stressMat;

                stressMat += dt * jaumannRate;
                m_stress[i][0] = stressMat(0, 0);
                m_stress[i][1] = stressMat(1, 1);
                m_stress[i][2] = stressMat(2, 2);
                m_stress[i][3] = stressMat(0, 1);
                m_stress[i][4] = stressMat(0, 2);
                m_stress[i][5] = stressMat(1, 2);
            }
        }
#ifndef ELASTICITY_DEBUG_SINGLE_THREAD
    }
#endif
}

void Elasticity_Test2026::computeForces()
{
    const unsigned int numParticles = m_model->numActiveParticles();
    if (numParticles == 0) return;

    Simulation *sim = Simulation::getCurrent();
    const unsigned int fluidModelIndex = m_model->getPointSetIndex();

#ifndef ELASTICITY_DEBUG_SINGLE_THREAD
    #pragma omp parallel default(shared)
    {
        #pragma omp for schedule(static)
#endif
        for (int i = 0; i < (int)numParticles; i++)
        {
            if (m_model->getParticleState(i) == ParticleState::Active)
            {
                const Vector3r &xi = m_model->getPosition(i);
                const Real rhoi = m_model->getDensity(i);
                const Real mi = m_model->getMass(i);

                // Use the deviatoric stress directly (no pressure subtraction!)
                const Vector6r &sigma_i = m_stress[i];

                Vector3r force = Vector3r::Zero();
                const unsigned int numNeighbors = sim->numberOfNeighbors(fluidModelIndex, fluidModelIndex, i);

                for (unsigned int j = 0; j < numNeighbors; j++)
                {
                    const unsigned int neighborIndex = sim->getNeighbor(fluidModelIndex, fluidModelIndex, i, j);
                    const Vector3r &xj = m_model->getPosition(neighborIndex);
                    const Real rhoj = m_model->getDensity(neighborIndex);
                    const Real mj = m_model->getMass(neighborIndex);

                    const Vector6r &sigma_j = m_stress[neighborIndex];

                    Vector3r gradW = sim->gradW(xi - xj);

                    // sigma_i/j * gradW_ij (matrix-vector product)
                    Vector3r sdi, sdj;
                    symMatTimesVec(sigma_i, gradW, sdi);
                    symMatTimesVec(sigma_j, gradW, sdj);

                    // Elastic force contribution (symmetric SPH formulation)
                    force += mi * mj * (sdi / (rhoi * rhoi) + sdj / (rhoj * rhoj));
                }

                // Add elastic acceleration to the particle
                Vector3r& ai = m_model->getAcceleration(i);
                ai += force / mi;
            }
        }
#ifndef ELASTICITY_DEBUG_SINGLE_THREAD
    }
#endif
}

void Elasticity_Test2026::reset()
{
    m_stress.resize(m_model->numActiveParticles(), Vector6r::Zero());
}

void Elasticity_Test2026::performNeighborhoodSearchSort()
{
}

void Elasticity_Test2026::saveState(BinaryFileWriter &binWriter)
{
    binWriter.writeBuffer((char*)m_stress.data(), m_stress.size() * sizeof(Vector6r));
}

void Elasticity_Test2026::loadState(BinaryFileReader &binReader)
{
    binReader.readBuffer((char*)m_stress.data(), m_stress.size() * sizeof(Vector6r));
}
