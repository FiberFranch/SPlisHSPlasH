#ifndef __Elasticity_Test2026_h__
#define __Elasticity_Test2026_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SPlisHSPlasH/NonPressureForceBase.h"

namespace SPH
{
	/** \brief This class implements a basic velocity-based SPH formulation for
	 * linearly elastic solids.
	 *
	 * References:
	 * - Yan et al. [2016]
	 * - Yang et al. [2017]
	 */
	class Elasticity_Test2026 : public NonPressureForceBase
	{
	protected:
		struct ElasticObject
		{
			std::vector<unsigned int> m_particleIndices;
			unsigned int m_nFixed = 0;

			ElasticObject() { m_nFixed = 0; }
			~ElasticObject() {}
		};
		
		Real m_youngsModulus;
		Real m_poissonRatio;
		Vector3r m_fixedBoxMin;
		Vector3r m_fixedBoxMax;
		std::vector<ElasticObject*> m_objects;

		// Deviatoric shear stress tensor for each particle
		std::vector<Vector6r> m_stress;

		// Parameters for tensile instability correction
		// TODO

		void computeStress();
		void computeForces();
		void determineFixedParticles();

		virtual void initParameters();
		/** This function is called after the simulation scene is loaded and all
		* parameters are initialized. While reading a scene file several parameters
		* can change. The deferred init function should initialize all values which
		* depend on these parameters.
		*/
		virtual void deferredInit();

		void initValues();
		void findObjects();

		//////////////////////////////////////////////////////////////////////////
		// multiplication of symmetric matrix, represented by a 6D vector, and a
		// 3D vector
		//////////////////////////////////////////////////////////////////////////
		FORCE_INLINE void symMatTimesVec(const Vector6r & M, const Vector3r & v, Vector3r &res)
		{
			res[0] = M[0] * v[0] + M[3] * v[1] + M[4] * v[2];
			res[1] = M[3] * v[0] + M[1] * v[1] + M[5] * v[2];
			res[2] = M[4] * v[0] + M[5] * v[1] + M[2] * v[2];
		}

	public:
		static std::string METHOD_NAME;
		static int YOUNGS_MODULUS;
		static int POISSON_RATIO;
		static int FIXED_BOX_MIN;
		static int FIXED_BOX_MAX;

		Elasticity_Test2026(FluidModel *model);
		virtual ~Elasticity_Test2026(void);

		static NonPressureForceBase* creator(FluidModel* model) { return new Elasticity_Test2026(model); }
		virtual std::string getMethodName() { return METHOD_NAME; }

		virtual void step();
		virtual void reset();
		virtual void performNeighborhoodSearchSort();

		virtual void saveState(BinaryFileWriter &binWriter);
		virtual void loadState(BinaryFileReader &binReader);
	};
}

#endif
