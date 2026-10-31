#include "Elasticity_Test2026.h"
#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "Utilities/Timing.h"

using namespace SPH;
using namespace GenParam;

std::string Elasticity_Test2026::METHOD_NAME = "Elasticity Test 2026";
int Elasticity_Test2026::YOUNGS_MODULUS = -1;
int Elasticity_Test2026::POISSON_RATIO = -1;

Elasticity_Test2026::Elasticity_Test2026(FluidModel *model) : NonPressureForceBase(model)
{
    const unsigned int numParticles = model->numActiveParticles();
    m_stress.resize(numParticles);

    m_youngsModulus = static_cast<Real>(100000.0);
    m_poissonRatio = static_cast<Real>(0.3);

    model->addField({ "deviatoric stress", METHOD_NAME, FieldType::Vector6, [&](const unsigned int i) -> Real* { return &m_stress[i][0]; }, true });
}

Elasticity_Test2026::~Elasticity_Test2026(void)
{
    m_model->removeFieldByName("deviatoric stress");
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

    #pragma omp parallel default(shared)
    {
        #pragma omp for schedule(static)
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
    }
}

void Elasticity_Test2026::computeForces()
{
    const unsigned int numParticles = m_model->numActiveParticles();
    if (numParticles == 0) return;

    Simulation *sim = Simulation::getCurrent();
    const unsigned int fluidModelIndex = m_model->getPointSetIndex();

    #pragma omp parallel default(shared)
    {
        #pragma omp for schedule(static)
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
    }
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
