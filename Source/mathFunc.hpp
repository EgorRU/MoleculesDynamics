#ifndef MATHFUNC_HPP
#define MATHFUNC_HPP

#include <QVector>
#include <QVector3D>
#include <cmath>

static inline QVector<QVector3D> computeLJForces(
	const QVector<QVector3D>& positions,
	double epsilon,
	double sigma,
	double weight
)
{
	int N = positions.size();
	QVector<QVector3D> accelerations(N, QVector3D(0, 0, 0));

	for (int i = 0; i < N; ++i)
	{
		for (int j = i + 1; j < N; ++j)
		{
			QVector3D rij = positions[i] - positions[j];
			double r = rij.length();
			if (r < 1e-6)
			{
				continue;
			}

			double sr = sigma / r;
			double sr6 = std::pow(sr, 6);
			double sr12 = sr6 * sr6;
			double F_over_m = (24.0 * epsilon / weight) * (2.0 * sr12 - sr6) / (r * r);
			QVector3D force = F_over_m * rij;

			accelerations[i] += force;
			accelerations[j] -= force;
		}
	}
	return accelerations;
}

// Метод Рунге–Кутта 4-го порядка (RK4):
// y_{n+1} = y_n + (dt / 6) * (k1 + 2k2 + 2k3 + k4)
// где:
// k1 = f(t_n, y_n)
// k2 = f(t_n + dt/2, y_n + dt/2 * k1)
// k3 = f(t_n + dt/2, y_n + dt/2 * k2)
// k4 = f(t_n + dt, y_n + dt * k3)
static inline void rungeKutta(
	QVector<QVector3D>& pos,
	QVector<QVector3D>& vel,
	double epsilon,
	double sigma,
	double weight,
	double dt
)
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

// Метод Верле (Verlet):
// x_{n+1} = 2x_n - x_{n-1} + a_n * dt^2
// v_n = (x_{n+1} - x_{n-1}) / (2 * dt)
static inline void verlet(
	QVector<QVector3D>& pos,
	QVector<QVector3D>& vel,
	QVector<QVector3D>& prev_pos,
	double epsilon,
	double sigma,
	double weight,
	double dt
)
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

// Метод скоростной Верле (Velocity Verlet):
// x_{n+1} = x_n + v_n * dt + (1/2) * a_n * dt^2
// a_{n+1} = f(x_{n+1})
// v_{n+1} = v_n + (1/2) * (a_n + a_{n+1}) * dt
static inline void velocityVerlet(
	QVector<QVector3D>& pos,
	QVector<QVector3D>& vel,
	double epsilon,
	double sigma,
	double weight,
	double dt
)
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

// Метод с перескоками (Leapfrog):
// v_{n+1/2} = v_n + (1/2) * a_n * dt
// x_{n+1} = x_n + v_{n+1/2} * dt
// v_{n+1} = v_{n+1/2} + (1/2) * a_{n+1} * dt
static inline void leapfrog(
	QVector<QVector3D>& pos,
	QVector<QVector3D>& vel,
	double epsilon,
	double sigma,
	double weight,
	double dt
)
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

// Метод Бимана-Шофилда (Beeman-Schofield):
// x_{n+1} = x_n + v_n * dt + (1/6) * (4a_n - a_{n-1}) * dt^2
// v_{n+1} = v_n + (1/6) * (2a_{n+1} + 5a_n - a_{n-1}) * dt
static inline void beemanSchofield(
	QVector<QVector3D>& pos,
	QVector<QVector3D>& vel,
	QVector<QVector3D>& prev_forces,
	double epsilon,
	double sigma,
	double weight,
	double dt
)
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

// Метод предиктор-корректор:
// Предиктор: x_{n+1}^p = x_n + v_n * dt + 0.5 * a_n * dt^2
//             v_{n+1}^p = v_n + a_n * dt
// Корректор: x_{n+1} = x_n + v_n * dt + 0.5 * a_{n+1}^p * dt^2
//            v_{n+1} = v_n + 0.5 * (a_n + a_{n+1}^p) * dt
static inline void predictorCorrector(
	QVector<QVector3D>& pos,
	QVector<QVector3D>& vel,
	QVector<QVector3D>& prev_forces,
	double epsilon,
	double sigma,
	double weight,
	double dt
)
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

static inline double reflect(double coord, double L, double& vel)
{
	double size = 2.0f * L;
	double shifted = coord + L;

	int n = static_cast<int>(shifted / size);
	double local = fmod(shifted, size);

	if (local < 0.0f)
	{
		local += size;
		n -= 1;
	}

	if (n % 2 != 0)
	{
		local = size - local;
		vel = -vel;
	}
	return local - L;
}

#endif