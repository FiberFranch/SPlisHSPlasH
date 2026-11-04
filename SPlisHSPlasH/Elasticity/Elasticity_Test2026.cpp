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
    LOG_INFO << "METHOD: " << METHOD_NAME;
    LOG_INFO << "youngsModulus: " << m_youngsModulus;
    LOG_INFO << "poissonsRatio: " << m_poissonRatio;
    LOG_INFO << "fixedBoxMin: " << m_fixedBoxMin.x() << ", " << m_fixedBoxMin.y() << ", " << m_fixedBoxMin.z();
    LOG_INFO << "fixedBoxMax: " << m_fixedBoxMax.x() << ", " << m_fixedBoxMax.y() << ", " << m_fixedBoxMax.z();
}

void Elasticity_Test2026::initParameters()
{
    NonPressureForceBase::initParameters();

    ParameterBase::GetFunc<Real> getFctYM = [&]()-> Real { return m_youngsModulus; };
    ParameterBase::SetFunc<Real> setFctYM = [&](Real val)
    {
        m_youngsModulus = val;
    };
    YOUNGS_MODULUS = createNumericParameter("youngsModulus", "Young`s modulus", getFctYM, setFctYM);
    setGroup(YOUNGS_MODULUS, "Fluid Model|Elasticity");
    setDescription(YOUNGS_MODULUS, "Stiffness of the elastic material");
    RealParameter* rparam = static_cast<RealParameter*>(getParameter(YOUNGS_MODULUS));
    rparam->setMinValue(0.0);

    ParameterBase::GetFunc<Real> getFctPR = [&]()-> Real { return m_poissonRatio; };
    ParameterBase::SetFunc<Real> setFctPR = [&](Real val)
    {
        m_poissonRatio = val;
    };
    POISSON_RATIO = createNumericParameter("poissonsRatio", "Poisson`s ratio", getFctPR, setFctPR);
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

    LOG_INFO << "Elasticity_Test2026::initParameters called";
}

/** Mark all particles in the bounding box as fixed.
*/
void Elasticity_Test2026::determineFixedParticles()
{
    const unsigned int numParticles = m_model->numActiveParticles();

    if (!m_fixedBoxMin.isZero() || !m_fixedBoxMax.isZero())
    {
        for (int i = 0; i < (int)numParticles; i++)
        {
            const Vector3r& x = m_model->getPosition0(i);
            if ((x[0] > m_fixedBoxMin[0]) && (x[1] > m_fixedBoxMin[1]) && (x[2] > m_fixedBoxMin[2]) &&
                (x[0] < m_fixedBoxMax[0]) && (x[1] < m_fixedBoxMax[1]) && (x[2] < m_fixedBoxMax[2]))
            {
                m_model->setParticleState(i, ParticleState::Fixed);
            }
        }
    }
}

void Elasticity_Test2026::initValues()
{
    const unsigned int numParticles = m_model->numActiveParticles();

    // Reset all particles to Active
#ifndef ELASTICITY_DEBUG_SINGLE_THREAD
    #pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < (int)numParticles; i++)
    {
        if (m_model->getParticleState(i) == ParticleState::Fixed)
            m_model->setParticleState(i, ParticleState::Active);
    }

    LOG_INFO << "Box min: " << m_fixedBoxMin.x() << ", " << m_fixedBoxMin.y() << ", " << m_fixedBoxMin.z();
    LOG_INFO << "Box max: " << m_fixedBoxMax.x() << ", " << m_fixedBoxMax.y() << ", " << m_fixedBoxMax.z();

    // mark all particles in the bounding box as fixed
    determineFixedParticles();
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

                Real maxStrain = strainRate.norm();
                Real maxSpin = spinRate.norm();
                LOG_INFO << "max |strainRate| = " << maxStrain << "  max |spinRate| = " << maxSpin;

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
        Real maxStress = 0.0;
        for (int i = 0; i < (int)numParticles; i++)
        {
            for (int k = 0; k < 6; k++)
                maxStress = std::max(maxStress, std::abs(m_stress[i][k]));
        }
        LOG_INFO << "Step " << TimeManager::getCurrent()->getTime() << "  max |stress| = " << maxStress;
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
