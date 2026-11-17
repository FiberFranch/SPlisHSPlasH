// #define ELASTICITY_DEBUG_SINGLE_THREAD

#include "Elasticity_Test2026.h"
#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/Utilities/MathFunctions.h"
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

    m_youngsModulus = static_cast<Real>(100000.0);
    m_poissonRatio = static_cast<Real>(0.3);
    m_fixedBoxMin.setZero();
    m_fixedBoxMax.setZero();

	m_restVolumes.resize(numParticles);
    m_current_to_initial_index.resize(numParticles);
	m_initial_to_current_index.resize(numParticles);

    m_L.resize(numParticles, Matrix3r::Zero());
    m_stress.resize(numParticles, Vector6r::Zero());

    model->addField({ "rest volume", METHOD_NAME, FieldType::Scalar, [&](const unsigned int i) -> Real* { return &m_restVolumes[i]; }, true });
	model->addField({ "correction matrix", METHOD_NAME, FieldType::Matrix3, [&](const unsigned int i) -> Real* { return &m_L[i](0,0); } });
    model->addField({ "deviatoric stress", METHOD_NAME, FieldType::Vector6, [&](const unsigned int i) -> Real* { return &m_stress[i][0]; }, true });

}

Elasticity_Test2026::~Elasticity_Test2026(void)
{
    m_model->removeFieldByName("rest volume");
	m_model->removeFieldByName("correction matrix");
    m_model->removeFieldByName("deviatoric stress");
}

void Elasticity_Test2026::deferredInit()
{
    initValues();
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

/** Initialize the particle neighborhoods in the reference configuration.
* Fix particles which lie in the user-defined bounding box.
* Finally, compute kernel gradient correction matrices and factorization.
*/
void Elasticity_Test2026::initValues()
{
    Simulation *sim = Simulation::getCurrent();
    NeighborhoodSearchWrapper* ns = sim->getNeighborhoodSearch();
    ns->findNeighbors();

    FluidModel *model = m_model;
    const unsigned int numParticles = model->numActiveParticles();
    const unsigned int fluidModelIndex = model->getPointSetIndex();

#ifndef ELASTICITY_DEBUG_SINGLE_THREAD
    #pragma omp parallel default(shared)
    {
        #pragma omp for schedule(static)
#endif
        for (int i = 0; i < (int)numParticles; i++)
        {
            m_current_to_initial_index[i] = i;
            m_initial_to_current_index[i] = i;

            // reset particle state
            if (m_model->getParticleState(i) == ParticleState::Fixed)
                m_model->setParticleState(i, ParticleState::Active);

            // compute rest volume
			Real density = model->getMass(i) * sim->W_zero();
            const Vector3r &xi0 = model->getPosition0(i);
            const unsigned int numNeighbors = sim->numberOfNeighbors(fluidModelIndex, fluidModelIndex, i);
            for (unsigned int j = 0; j < numNeighbors; j++)
            {
                const unsigned int neighborIndex0 = sim->getNeighbor(fluidModelIndex, fluidModelIndex, i, j);
                const Vector3r &xj0 = model->getPosition0(neighborIndex0);
                density += model->getMass(neighborIndex0) * sim->W(xi0 - xj0);
            }
            m_restVolumes[i] = model->getMass(i) / density;
            m_stress[i].setZero();
        }
#ifndef ELASTICITY_DEBUG_SINGLE_THREAD
    }
#endif

    // mark all particles in the bounding box as fixed
    determineFixedParticles();

    // compute kernel gradient correction matrix
	computeMatrixL();
}

void Elasticity_Test2026::step()
{
    computeMatrixL();
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

                    // Compute correction matrix Li from Yan paper
                    gradV += Vj * (vj - vi) * (m_L[i] * sim->gradW(xi - xj)).transpose();
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

                // Compute Jaumann stress rate from Yan eq 12
                Matrix3r jaumannRate = spinRate * stressMat - stressMat * spinRate + static_cast<Real>(2.0) * mu * strainRate;

                stressMat += dt * jaumannRate;
                m_stress[i][0] = stressMat(0, 0);
                m_stress[i][1] = stressMat(1, 1);
                m_stress[i][2] = stressMat(2, 2);
                m_stress[i][3] = stressMat(0, 1);
                m_stress[i][4] = stressMat(0, 2);
                m_stress[i][5] = stressMat(1, 2);

                if ((i == 0 || i == 4071) && TimeManager::getCurrent()->getTime() < 0.01)
                {
#ifndef ELASTICITY_DEBUG_SINGLE_THREAD
                    #pragma omp critical(logElasticityDiag)
#endif
                    {
                        LOG_INFO << "############ STEP t=" << TimeManager::getCurrent()->getTime()
                                << "  i=" << i << " ############";

                        LOG_INFO << "xi        = " << xi.transpose();
                        LOG_INFO << "vi        = " << vi.transpose();
                        LOG_INFO << "gradV row0= " << gradV.row(0);
                        LOG_INFO << "gradV row1= " << gradV.row(1);
                        LOG_INFO << "gradV row2= " << gradV.row(2);
                        LOG_INFO << "strainRate row0= " << strainRate.row(0);
                        LOG_INFO << "strainRate row1= " << strainRate.row(1);
                        LOG_INFO << "strainRate row2= " << strainRate.row(2);
                        LOG_INFO << "spinRate   row0= " << spinRate.row(0);
                        LOG_INFO << "spinRate   row1= " << spinRate.row(1);
                        LOG_INFO << "spinRate   row2= " << spinRate.row(2);
                        LOG_INFO << "stress     row0= " << stressMat.row(0);
                        LOG_INFO << "stress     row1= " << stressMat.row(1);
                        LOG_INFO << "stress     row2= " << stressMat.row(2);

                        if (numNeighbors > 0)
                        {
                            const unsigned int nb = sim->getNeighbor(fluidModelIndex, fluidModelIndex, i, 0);
                            const Vector3r &xj = m_model->getPosition(nb);
                            const Vector3r &vj = m_model->getVelocity(nb);
                            const Real Vj = m_model->getMass(nb) / m_model->getDensity(nb);
                            const Vector3r xij = xi - xj;
                            const Vector3r vij = vi - vj;

                            LOG_INFO << "neighbor j= " << nb;
                            LOG_INFO << "xj        = " << xj.transpose();
                            LOG_INFO << "vj        = " << vj.transpose();
                            LOG_INFO << "xi - xj   = " << xij.transpose();
                            LOG_INFO << "vi - vj   = " << vij.transpose();
                            LOG_INFO << "Vj        = " << Vj;
                            LOG_INFO << "gradW     = " << sim->gradW(xij).transpose();
                        }

                        LOG_INFO << "############ END i=" << i << " ############";
#ifndef ELASTICITY_DEBUG_SINGLE_THREAD
                    }
#endif
                }
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
    initValues();
}

void Elasticity_Test2026::performNeighborhoodSearchSort()
{
    const unsigned int numPart = m_model->numActiveParticles();
	if (numPart == 0)
		return;

    Simulation *sim = Simulation::getCurrent();
	NeighborhoodSearchWrapper* ns = sim->getNeighborhoodSearch();
	ns->applyZSort(m_model->getPointSetIndex(), &m_current_to_initial_index[0]);
	ns->applyZSort(m_model->getPointSetIndex(), &m_L[0]);
	ns->applyZSort(m_model->getPointSetIndex(), &m_restVolumes[0]);
	ns->applyZSort(m_model->getPointSetIndex(), &m_stress[0]);

    for (unsigned int i = 0; i < numPart; i++)
		m_initial_to_current_index[m_current_to_initial_index[i]] = i;
}

void Elasticity_Test2026::computeMatrixL()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int numParticles = m_model->numActiveParticles();
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();

#ifndef ELASTICITY_DEBUG_SINGLE_THREAD
    #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static)
#endif
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r &xi = m_model->getPosition(i);
			Matrix3r L;
			L.setZero();

			const unsigned int numNeighbors = sim->numberOfNeighbors(fluidModelIndex, fluidModelIndex, i);

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < numNeighbors; j++)
			{
				const unsigned int neighborIndex = sim->getNeighbor(fluidModelIndex, fluidModelIndex, i, j);

				const Vector3r &xj = m_model->getPosition(neighborIndex);
				const Vector3r xj_xi = xj - xi;
				const Vector3r gradW = sim->gradW(xj_xi);
				const Real Vj = m_model->getMass(neighborIndex) / m_model->getDensity(neighborIndex);

				// minus because gradW(xij) == -gradW(xji)
				L -= Vj * gradW * xj_xi.transpose();
			}

			// add 1 to z-component. otherwise we get a singular matrix in 2D
			if (sim->is2DSimulation())
				L(2, 2) = 1.0;

			bool invertible = false;
			L.computeInverseWithCheck(m_L[i], invertible, static_cast<Real>(1e-9));
			if (!invertible)
			{
				LOG_INFO << "Matrix not invertible.";
				MathFunctions::pseudoInverse(L, m_L[i]);
				//m_L[i] = Matrix3r::Identity();
			}
		}
#ifndef ELASTICITY_DEBUG_SINGLE_THREAD
    }
#endif
}

void Elasticity_Test2026::saveState(BinaryFileWriter &binWriter)
{
    binWriter.writeBuffer((char*)m_stress.data(), m_stress.size() * sizeof(Vector6r));
}

void Elasticity_Test2026::loadState(BinaryFileReader &binReader)
{
    binReader.readBuffer((char*)m_stress.data(), m_stress.size() * sizeof(Vector6r));
}
