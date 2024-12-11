#pragma once

#include <core/array.h>

M_NAMESPACE_BEGIN
	/**
	 * \brief Calculates the Fresnel reflection coefficient
	 * at a planar interface between two dielectrics
	 *
	 * \param cos_theta_i
	 *      Cosine of the angle between the surface normal and the incident ray
	 *
	 * \param eta
	 *      Relative refractive index of the interface. A value greater than 1.0
	 *      means that the surface normal is pointing into the region of lower
	 *      density.
	 *
	 * \return A tuple (F, cos_theta_t, eta_it, eta_ti) consisting of
	 *
	 *     F           Fresnel reflection coefficient.
	 *
	 *     cos_theta_t Cosine of the angle between the surface normal and the transmitted ray
	 *
	 *     eta_it      Relative index of refraction in the direction of travel.
	 *
	 *     eta_ti      Reciprocal of the relative index of refraction in the
	 *                 direction of travel. This also happens to be equal to
	 *                 the scale factor that must be applied to the X and Y
	 *                 component of the refracted direction.
	 */
	template <typename Scalar>
	std::tuple<Scalar, Scalar, Scalar, Scalar> fresnel(Scalar cos_theta_i, Scalar eta) {
		auto outside_mask = cos_theta_i >= 0.f;

		Scalar rcp_eta = Scalar(1) / eta,
		       eta_it = outside_mask ? eta : rcp_eta,
		       eta_ti = outside_mask ? rcp_eta : eta;

		/* Using Snell's law, calculate the squared sine of the
		   angle between the surface normal and the transmitted ray */
		Scalar cos_theta_t_sqr = Scalar(1) - eta_ti * eta_ti * (Scalar(1) - cos_theta_i * cos_theta_i);

		/* Find the absolute cosines of the incident/transmitted rays */
		Scalar cos_theta_i_abs = abs(cos_theta_i);
		Scalar cos_theta_t_abs = safe_sqrt(cos_theta_t_sqr);

		auto index_matched = eta == Scalar(1);
		auto special_case = index_matched || cos_theta_i_abs == Scalar(0);

		Scalar r_sc = index_matched ? Scalar(0) : Scalar(1);

		/* Amplitudes of reflected waves */
		Scalar a_s = (cos_theta_i_abs - eta_it * cos_theta_t_abs) / (eta_it * cos_theta_t_abs + cos_theta_i_abs);
		Scalar a_p = (cos_theta_t_abs - eta_it * cos_theta_i_abs) / (eta_it * cos_theta_i_abs + cos_theta_t_abs);

		Scalar r = Scalar(0.5) * (a_s * a_s + a_p * a_p);
		if (special_case) {
			r = r_sc;
		}

		/* Adjust the sign of the transmitted direction */
		Scalar cos_theta_t = -mulsign(cos_theta_t_abs, cos_theta_i);

		return {r, cos_theta_t, eta_it, eta_ti};
	}

	/**
	 * \brief Calculates the Fresnel reflection coefficient at a planar
	 * interface of a conductor, i.e. a surface with a complex-valued relative index
	 * of refraction
	 *
	 * \remark
	 *      The implementation assumes that cos_theta_i > 0, i.e. light enters
	 *      from *outside* of the conducting layer (generally a reasonable
	 *      assumption unless very thin layers are being simulated)
	 *
	 * \param cos_theta_i
	 *      Cosine of the angle between the surface normal and the incident ray
	 *
	 * \param eta
	 *      Relative refractive index (complex-valued)
	 *
	 * \return The Fresnel reflection coefficient.
	 */
	template <typename Scalar>
	Scalar fresnel_conductor(Scalar cos_theta_i, Scalar eta, Scalar k) {
		Scalar cos_theta_i_2 = cos_theta_i * cos_theta_i,
		       sin_theta_i_2 = Scalar(1) - cos_theta_i_2,
		       sin_theta_i_4 = sin_theta_i_2 * sin_theta_i_2;

		auto eta_r = eta;
		auto eta_i = k;

		Scalar temp_1 = eta_r * eta_r - eta_i * eta_i - sin_theta_i_2,
		       a_2_pb_2 = safe_sqrt(temp_1 * temp_1 + Scalar(4) * eta_i * eta_i * eta_r * eta_r),
		       a = safe_sqrt(Scalar(0.5) * (a_2_pb_2 + temp_1));

		Scalar term_1 = a_2_pb_2 + cos_theta_i_2;
		Scalar term_2 = Scalar(2) * cos_theta_i * a;

		Scalar r_s = (term_1 - term_2) / (term_1 + term_2);

		Scalar term_3 = a_2_pb_2 * cos_theta_i_2 + sin_theta_i_4;
		Scalar term_4 = term_2 * sin_theta_i_2;

		Scalar r_p = r_s * (term_3 - term_4) / (term_3 + term_4);

		return Scalar(0.5) * (r_s + r_p);
	}

	// Reflection in local coordinates
	template <typename Scalar, DeviceType Device>
	TArray<Scalar, 3, ArrayType::Vector, Device> reflect(const TArray<Scalar, 3, ArrayType::Vector, Device>& wi) {
		return TArray<Scalar, 3, ArrayType::Vector, Device>({-wi.x(), -wi.y(), wi.z()});
	}

	// Reflect \c wi with respect to a given surface normal
	template <typename Scalar, DeviceType Device>
	TArray<Scalar, 3, ArrayType::Vector, Device> reflect(const TArray<Scalar, 3, ArrayType::Vector, Device>& wi,
	                                                     const TArray<Scalar, 3, ArrayType::Normal, Device>& m) {
		return TArray<Scalar, 3, ArrayType::Vector, Device>(m) * Scalar(2) * wi.dot(m) - wi;
	}

	/**
	 * \brief Refraction in local coordinates
	 *
	 * The 'cos_theta_t' and 'eta_ti' parameters are given by the last two tuple
	 * entries returned by the \ref fresnel and \ref fresnel_polarized functions.
	 */
	template <typename Scalar, DeviceType Device>
	TArray<Scalar, 3, ArrayType::Vector, Device> refract(const TArray<Scalar, 3, ArrayType::Vector, Device>& wi,
	                                                     Scalar cos_theta_t, Scalar eta_ti) {
		return TArray<Scalar, 3, ArrayType::Vector, Device>({-eta_ti * wi.x(), -eta_ti * wi.y(), cos_theta_t});
	}

	/**
	 * \brief Refract \c wi with respect to a given surface normal
	 *
	 * \param wi
	 *     Direction to refract
	 * \param m
	 *     Surface normal
	 * \param cos_theta_t
	 *     Cosine of the angle between the normal the transmitted
	 *     ray, as computed e.g. by \ref fresnel.
	 * \param eta_ti
	 *     Relative index of refraction (transmitted / incident)
	 */
	template <typename Scalar, DeviceType Device>
	TArray<Scalar, 3, ArrayType::Vector, Device> refract(const TArray<Scalar, 3, ArrayType::Vector, Device>& wi,
	                                                     const TArray<Scalar, 3, ArrayType::Normal, Device>& m,
	                                                     Scalar cos_theta_t, Scalar eta_ti) {
		return m * (wi.dot(m) * eta_ti + cos_theta_t) - wi * eta_ti;
	}

M_NAMESPACE_END
