#include "MoleculesDynamics.h"

QString labels[6] =
{
	"Рунге-Кутта 4 порядка",
	"Верле",
	"Скоростной Верле",
	"С перескоками",
	"Бимана-Шофилда",
	"Предиктор-корректор"
};

MoleculesDynamics::MoleculesDynamics(QWidget* parent)
	: QMainWindow(parent), currentStep(0)
{
	animationTimer = new QTimer();
	gridLayout = new QGridLayout();
	currentTimerLabel = new QLabel();

	NSpinBox = new QDoubleSpinBox();
	epsilonSpinBox = new QDoubleSpinBox();
	sigmaSpinBox = new QDoubleSpinBox();
	weightSpinBox = new QDoubleSpinBox();
	densitySpinBox = new QDoubleSpinBox();
	stepSpinBox = new QDoubleSpinBox();
	speedSpinBox = new QDoubleSpinBox();

	MSEComboBox = new QComboBox();

	rungeKuttaLabel = new QLabel();
	verletLabel = new QLabel();
	velocityVerletLabel = new QLabel();
	leapfrogLabel = new QLabel();
	beemanSchofieldLabel = new QLabel();
	predictorCorrectorLabel = new QLabel();

	spinBoxConfigs =
	{
		{stepSpinBox, {100, 100000, 100, 10000, "Шагов моделирования"}},
		{NSpinBox, {2, 1000, 1, 10, "Количество молекул"}},
		{densitySpinBox, {0.1, 5, 0.1, 0.8, "Плотность (ρ)"}},
		{weightSpinBox, {0.1, 1000, 1, 1, "Масса"}},
		{epsilonSpinBox, {0.1, 1000, 1, 1, "ε"}},
		{sigmaSpinBox, {0.1, 1000, 1, 1, "σ"}},
		{speedSpinBox, {0.01, 3, 0.01, 0.2, "Диапазон начальных скорост"}},
	};

	labelsMSE = new QLabel * [6]
	{
		rungeKuttaLabel,
		verletLabel,
		velocityVerletLabel,
		leapfrogLabel,
		beemanSchofieldLabel,
		predictorCorrectorLabel,
	};

	setupUI();
}

void MoleculesDynamics::setupUI()
{
	centralWidget = new QWidget(this);
	setCentralWidget(centralWidget);

	animationTimer->setInterval(10);
	connect(animationTimer, &QTimer::timeout, this, &MoleculesDynamics::animateScatters);

	QHBoxLayout* mainLayout = new QHBoxLayout(centralWidget);
	centralWidget->setLayout(mainLayout);

	QWidget* plotsWidget = new QWidget(this);
	QVBoxLayout* plotsLayout = new QVBoxLayout(plotsWidget);
	plotsLayout->addLayout(gridLayout);
	mainLayout->addWidget(plotsWidget, 1);

	QWidget* controlsWidget = new QWidget(this);
	QVBoxLayout* controlsLayout = new QVBoxLayout(controlsWidget);
	controlsLayout->setAlignment(Qt::AlignTop | Qt::AlignLeft);

	QHBoxLayout* currentStepLayout = new QHBoxLayout();
	currentStepLayout->setAlignment(Qt::AlignLeft);
	currentStepLayout->addWidget(currentTimerLabel);
	controlsLayout->addLayout(currentStepLayout);

	for (auto it = spinBoxConfigs.begin(); it != spinBoxConfigs.end(); ++it) {
		QDoubleSpinBox* spinBox = it.key();
		const SpinBoxConfig& config = it.value();

		spinBox->setMinimum(config.min);
		spinBox->setMaximum(config.max);
		spinBox->setSingleStep(config.step);
		spinBox->setValue(config.value);
		spinBox->setFixedWidth(80);

		QHBoxLayout* layout = new QHBoxLayout();
		layout->setAlignment(Qt::AlignLeft);
		QLabel* label = new QLabel(config.label);
		label->setFixedWidth(180);
		layout->addWidget(label);
		layout->addWidget(spinBox);
		layout->addStretch();

		controlsLayout->addLayout(layout);

		connect(spinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
			this, &MoleculesDynamics::resetAnimation);
	}

	QHBoxLayout* MSELayout = new QHBoxLayout();
	MSELayout->setAlignment(Qt::AlignLeft);
	QLabel* MSELabel = new QLabel("MSE относительно");
	MSELabel->setFixedWidth(180);
	for (int i = 0; i < 6; ++i)
	{
		MSEComboBox->addItem(labels[i]);
	}
	MSELayout->addWidget(MSELabel);
	MSELayout->addWidget(MSEComboBox);
	MSELayout->addStretch();
	controlsLayout->addLayout(MSELayout);
	connect(MSEComboBox, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &MoleculesDynamics::resetAnimation);

	controlsLayout->addWidget(rungeKuttaLabel);
	controlsLayout->addWidget(verletLabel);
	controlsLayout->addWidget(velocityVerletLabel);
	controlsLayout->addWidget(leapfrogLabel);
	controlsLayout->addWidget(beemanSchofieldLabel);
	controlsLayout->addWidget(predictorCorrectorLabel);

	QPushButton* restartButton = new QPushButton("Начать заново");
	restartButton->setFixedWidth(180);
	controlsLayout->addWidget(restartButton);
	connect(restartButton, &QPushButton::clicked, this, &MoleculesDynamics::resetAnimation);

	mainLayout->addWidget(controlsWidget, 0);

	for (int i = 0; i < 6; ++i)
	{
		createScatter(i);
		int row = i / 3;
		int col = i % 3;

		QWidget* plotContainer = new QWidget();
		QVBoxLayout* plotLayout = new QVBoxLayout(plotContainer);
		plotLayout->setContentsMargins(0, 0, 0, 0);
		plotLayout->setSpacing(2);

		QLabel* plotLabel = new QLabel(labels[i]);
		plotLabel->setAlignment(Qt::AlignCenter);
		plotLayout->addWidget(plotLabel);
		plotLayout->addWidget(scatterContainers[i]);

		gridLayout->addWidget(plotContainer, row, col);
	}
	resetAnimation();
}

void MoleculesDynamics::createScatter(int index)
{
	scatters[index] = new Q3DScatter();
	scatterContainers[index] = QWidget::createWindowContainer(scatters[index], centralWidget);
	scatterContainers[index]->setMinimumSize(300, 300);
	scatterContainers[index]->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
	scatterDataProxies[index] = new QScatterDataProxy();
	scatterSeries[index] = new QScatter3DSeries(scatterDataProxies[index]);
	scatterSeries[index]->setItemSize(0.1f);
	scatterSeries[index]->setBaseColor(Qt::blue);
	scatters[index]->addSeries(scatterSeries[index]);
	scatters[index]->axisX()->setTitle("X");
	scatters[index]->axisY()->setTitle("Y");
	scatters[index]->axisZ()->setTitle("Z");
	double N = NSpinBox->value();
	double L = pow(N / densitySpinBox->value(), 1.0 / 3.0);
	scatters[index]->axisX()->setRange(-L, L);
	scatters[index]->axisY()->setRange(-L, L);
	scatters[index]->axisZ()->setRange(-L, L);
	scatters[index]->setShadowQuality(QAbstract3DGraph::ShadowQualityNone);
	scatters[index]->scene()->activeCamera()->setCameraPreset(Q3DCamera::CameraPresetFront);
	scatters[index]->scene()->activeCamera()->setZoomLevel(120);
}

void MoleculesDynamics::animateScatters()
{
	int totalSteps = stepSpinBox->value();
	if (currentStep >= totalSteps)
	{
		resetAnimation();
		return;
	}

	double epsilon = epsilonSpinBox->value();
	double sigma = sigmaSpinBox->value();
	double weight = weightSpinBox->value();

	for (int i = 0; i < 6; ++i)
	{
		switch (i)
		{
		case 0:
			rungeKutta(method_positions[i], method_velocities[i], epsilon, sigma, weight);
			break;
		case 1:
			verlet(method_positions[i], method_velocities[i], verlet_prev_positions, epsilon, sigma, weight);
			break;
		case 2:
			velocityVerlet(method_positions[i], method_velocities[i], epsilon, sigma, weight);
			break;
		case 3:
			leapfrog(method_positions[i], method_velocities[i], epsilon, sigma, weight);
			break;
		case 4:
			beemanSchofield(method_positions[i], method_velocities[i], method_prev_forces[i], epsilon, sigma, weight);
			break;
		case 5:
			predictorCorrector(method_positions[i], method_velocities[i], method_prev_forces[i], epsilon, sigma, weight);
			break;
		}
		double L = pow(NSpinBox->value() / densitySpinBox->value(), 1.0 / 3.0);
		applyBoundaryConditions(method_positions[i], method_velocities[i], L);
	}

	currentStep++;
	int referenceMethod = MSEComboBox->currentIndex();
	for (int i = 0; i < 6; ++i)
	{
		double currentMSE = calculateMSE(method_positions[i], method_positions[referenceMethod]);
		accumulatedMSE[i] += currentMSE;
	}

	updateScatters();
	updateMSELabels();
	currentTimerLabel->setText(QString("Шаг моделирования: %1 из %2").arg(currentStep).arg(totalSteps));
}

void MoleculesDynamics::updateScatters()
{
	int N = NSpinBox->value();
	double L = pow(N / densitySpinBox->value(), 1.0 / 3.0);

	for (int i = 0; i < 6; ++i)
	{
		QScatterDataArray* dataArray = new QScatterDataArray();
		dataArray->resize(N);
		for (int j = 0; j < N; ++j)
		{
			(*dataArray)[j].setPosition(method_positions[i][j]);
		}
		scatterDataProxies[i]->resetArray(dataArray);
	}

	for (int i = 0; i < 6; ++i)
	{
		scatters[i]->axisX()->setRange(-L, L);
		scatters[i]->axisY()->setRange(-L, L);
		scatters[i]->axisZ()->setRange(-L, L);
	}
}

void MoleculesDynamics::resetAnimation()
{
	currentStep = 0;
	accumulatedMSE.clear();
	accumulatedMSE.resize(6);
	for (double& mse : accumulatedMSE)
	{
		mse = 0.0;
	}

	int N = NSpinBox->value();
	double L = pow(N / densitySpinBox->value(), 1.0 / 3.0);
	double speed = speedSpinBox->value();

	positions.clear();
	velocities.clear();
	positions.resize(N);
	velocities.resize(N);

	for (int i = 0; i < N; ++i)
	{
		double x = -L + QRandomGenerator::global()->generateDouble() * 2 * L;
		double y = -L + QRandomGenerator::global()->generateDouble() * 2 * L;
		double z = -L + QRandomGenerator::global()->generateDouble() * 2 * L;
		positions[i] = QVector3D(x, y, z);
	}

	for (int i = 0; i < N; ++i)
	{
		double vx = -speed + QRandomGenerator::global()->generateDouble() * speed * 2;
		double vy = -speed + QRandomGenerator::global()->generateDouble() * speed * 2;
		double vz = -speed + QRandomGenerator::global()->generateDouble() * speed * 2;
		velocities[i] = QVector3D(vx, vy, vz);
	}

	QVector3D totalMomentum(0, 0, 0);
	for (const auto& v : velocities)
	{
		totalMomentum += v;
	}
	totalMomentum /= N;
	for (auto& v : velocities)
	{
		v -= totalMomentum;
	}

	method_positions.clear();
	method_velocities.clear();
	method_prev_forces.clear();

	for (int i = 0; i < 6; ++i)
	{
		method_positions.push_back(positions);
		method_velocities.push_back(velocities);
		method_prev_forces.push_back(QVector<QVector3D>(N, QVector3D(0, 0, 0)));
	}

	verlet_prev_positions = method_positions[1];
	for (int j = 0; j < N; ++j)
	{
		verlet_prev_positions[j] -= method_velocities[1][j] * dt;
	}

	double epsilon = epsilonSpinBox->value();
	double sigma = sigmaSpinBox->value();
	double weight = weightSpinBox->value();

	for (int i = 4; i <= 5; ++i)
	{
		method_prev_forces[i] = computeLJForces(method_positions[i], epsilon, sigma, weight);
	}

	updateScatters();
	updateMSELabels();
	animationTimer->start();
}

void MoleculesDynamics::applyBoundaryConditions(QVector<QVector3D>& pos, QVector<QVector3D>& vel, double L)
{
	double boxSize = 2.0 * L;
	for (int i = 0; i < pos.size(); ++i)
	{
		double x = pos[i].x();
		double y = pos[i].y();
		double z = pos[i].z();

		x = -L + fmod((x + L), boxSize);
		y = -L + fmod((y + L), boxSize);
		z = -L + fmod((z + L), boxSize);

		pos[i].setX(x);
		pos[i].setY(y);
		pos[i].setZ(z);
	}
}

double MoleculesDynamics::calculateMSE(const QVector<QVector3D>& pos1, const QVector<QVector3D>& pos2)
{
	double mse = 0.0;
	for (int i = 0; i < pos1.size(); ++i)
	{
		QVector3D diff = pos1[i] - pos2[i];
		mse += diff.lengthSquared();
	}
	return mse / pos1.size();
}

void MoleculesDynamics::updateMSELabels()
{
	int referenceMethod = MSEComboBox->currentIndex();
	for (int i = 0; i < 6; ++i)
	{
		double averageMSE = accumulatedMSE[i] / currentStep;
		QString mseText = QString("MSE %1: %2").arg(labels[i]).arg(averageMSE, 0, 'f', 6);
		labelsMSE[i]->setText(mseText);
	}
}