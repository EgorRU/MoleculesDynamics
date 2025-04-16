#ifndef MATHFUNC_HPP
#define MATHFUNC_HPP

#include <QVector>
#include <QVector3D>
#include <cmath>

static const double dt = 0.00001;

static inline QVector<QVector3D> computeLJForces(const QVector<QVector3D>& positions, double epsilon, double sigma, double weight)
{
	int N = positions.size();
	QVector<QVector3D> accelerations(N, QVector3D(0, 0, 0));

	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			if (i == j) continue;

			QVector3D rij = positions[i] - positions[j];
			double r = rij.length();
			if (r < 1e-6) continue;

			double sr = sigma / r;
			double sr6 = std::pow(sr, 6);
			double sr12 = sr6 * sr6;
			double F_over_m = (24.0 * epsilon / weight) * (2.0 * sr12 - sr6) / (r * r);
			accelerations[i] += F_over_m * rij;
		}
	}
	return accelerations;
}

static inline void rungeKutta(QVector<QVector3D>& pos, QVector<QVector3D>& vel, double epsilon, double sigma, double weight)
{
	int N = pos.size();
	QVector<QVector3D> k1_x(N), k1_v(N), k2_x(N), k2_v(N), k3_x(N), k3_v(N), k4_x(N), k4_v(N);
	QVector<QVector3D> temp_pos(N), temp_vel(N);

	k1_x = vel;
	k1_v = computeLJForces(pos, epsilon, sigma, weight);

	for (int i = 0; i < N; ++i)
	{
		temp_pos[i] = pos[i] + 0.5 * dt * k1_x[i];
		temp_vel[i] = vel[i] + 0.5 * dt * k1_v[i];
	}
	k2_x = temp_vel;
	k2_v = computeLJForces(temp_pos, epsilon, sigma, weight);

	for (int i = 0; i < N; ++i)
	{
		temp_pos[i] = pos[i] + 0.5 * dt * k2_x[i];
		temp_vel[i] = vel[i] + 0.5 * dt * k2_v[i];
	}
	k3_x = temp_vel;
	k3_v = computeLJForces(temp_pos, epsilon, sigma, weight);

	for (int i = 0; i < N; ++i)
	{
		temp_pos[i] = pos[i] + dt * k3_x[i];
		temp_vel[i] = vel[i] + dt * k3_v[i];
	}
	k4_x = temp_vel;
	k4_v = computeLJForces(temp_pos, epsilon, sigma, weight);

	for (int i = 0; i < N; ++i)
	{
		pos[i] += (dt / 6.0) * (k1_x[i] + 2.0 * k2_x[i] + 2.0 * k3_x[i] + k4_x[i]);
		vel[i] += (dt / 6.0) * (k1_v[i] + 2.0 * k2_v[i] + 2.0 * k3_v[i] + k4_v[i]);
	}
}

static inline void verlet(QVector<QVector3D>& pos, QVector<QVector3D>& vel, QVector<QVector3D>& prev_pos, double epsilon, double sigma, double weight)
{
	int N = pos.size();
	QVector<QVector3D> forces = computeLJForces(pos, epsilon, sigma, weight);

	for (int i = 0; i < N; ++i)
	{
		QVector3D r_new = 2.0 * pos[i] - prev_pos[i] + forces[i] * dt * dt;
		vel[i] = (r_new - prev_pos[i]) / (2.0 * dt);
		prev_pos[i] = pos[i];
		pos[i] = r_new;
	}
}

static inline void velocityVerlet(QVector<QVector3D>& pos, QVector<QVector3D>& vel, double epsilon, double sigma, double weight)
{
	int N = pos.size();
	QVector<QVector3D> forces = computeLJForces(pos, epsilon, sigma, weight);

	for (int i = 0; i < N; ++i)
	{
		pos[i] += vel[i] * dt + 0.5 * forces[i] * dt * dt;
	}

	QVector<QVector3D> new_forces = computeLJForces(pos, epsilon, sigma, weight);

	for (int i = 0; i < N; ++i)
	{
		vel[i] += 0.5 * (forces[i] + new_forces[i]) * dt;
	}
}

static inline void leapfrog(QVector<QVector3D>& pos, QVector<QVector3D>& vel, double epsilon, double sigma, double weight)
{
	int N = pos.size();
	QVector<QVector3D> forces = computeLJForces(pos, epsilon, sigma, weight);
	QVector<QVector3D> half_vel(N);

	for (int i = 0; i < N; ++i)
	{
		half_vel[i] = vel[i] + 0.5 * forces[i] * dt;
	}

	for (int i = 0; i < N; ++i)
	{
		pos[i] += half_vel[i] * dt;
	}

	forces = computeLJForces(pos, epsilon, sigma, weight);
	for (int i = 0; i < N; ++i)
	{
		vel[i] = half_vel[i] + 0.5 * forces[i] * dt;
	}
}

static inline void beemanSchofield(QVector<QVector3D>& pos, QVector<QVector3D>& vel, QVector<QVector3D>& prev_forces, double epsilon, double sigma, double weight)
{
	int N = pos.size();
	QVector<QVector3D> forces = computeLJForces(pos, epsilon, sigma, weight);

	for (int i = 0; i < N; ++i)
	{
		pos[i] += vel[i] * dt + (dt * dt / 6.0) * (4.0 * forces[i] - prev_forces[i]);
	}

	QVector<QVector3D> new_forces = computeLJForces(pos, epsilon, sigma, weight);

	for (int i = 0; i < N; ++i)
	{
		vel[i] += (dt / 6.0) * (2.0 * new_forces[i] + 5.0 * forces[i] - prev_forces[i]);
	}

	prev_forces = forces;
}

static inline void predictorCorrector(QVector<QVector3D>& pos, QVector<QVector3D>& vel, QVector<QVector3D>& prev_forces, double epsilon, double sigma, double weight)
{
	int N = pos.size();
	QVector<QVector3D> predicted_pos(N), predicted_vel(N);

	for (int i = 0; i < N; ++i)
	{
		predicted_pos[i] = pos[i] + vel[i] * dt + 0.5 * prev_forces[i] * dt * dt;
		predicted_vel[i] = vel[i] + prev_forces[i] * dt;
	}

	QVector<QVector3D> predicted_forces = computeLJForces(predicted_pos, epsilon, sigma, weight);

	for (int i = 0; i < N; ++i)
	{
		pos[i] += vel[i] * dt + 0.5 * predicted_forces[i] * dt * dt;
		vel[i] += 0.5 * (prev_forces[i] + predicted_forces[i]) * dt;
	}

	prev_forces = predicted_forces;
}

static inline float reflect(float coord, float L, float& vel)
{
	float size = 2.0f * L;
	float shifted = coord + L;

	int n = static_cast<int>(shifted / size);
	float local = shifted - n * size;

	if (local < 0.0f)
	{
		local += size;
		n -= 1;
	}

	bool isEven = (n % 2 == 0);
	if (!isEven)
	{
		local = size - local;
		vel = -vel;
	}
	return local - L;
}

#endif