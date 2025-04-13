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
	: QMainWindow(parent), currentStep(0), animationSpeed(1.0)
{
	animationTimer = new QTimer();
	gridLayout = new QGridLayout();
	currentTimerLabel = new QLabel();
	NSpinBox = new QDoubleSpinBox();
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

	setupUI();
}

void MoleculesDynamics::setupUI()
{
	QMenuBar* menuBar = new QMenuBar(this);
	setMenuBar(menuBar);

	QMenu* viewMenu = menuBar->addMenu("Внешний вид");

	QAction* resetCameraAction = new QAction("Сбросить камеры", this);
	viewMenu->addAction(resetCameraAction);
	connect(resetCameraAction, &QAction::triggered, this, &MoleculesDynamics::resetCamera);

	centralWidget = new QWidget(this);
	setCentralWidget(centralWidget);

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

	QHBoxLayout* NLayout = new QHBoxLayout();
	NLayout->setAlignment(Qt::AlignLeft);
	QLabel* NLabel = new QLabel("Количество молекул");
	NLabel->setFixedWidth(180);
	NSpinBox->setMinimum(2);
	NSpinBox->setMaximum(1000);
	NSpinBox->setValue(10);
	NSpinBox->setFixedWidth(80);
	NLayout->addWidget(NLabel);
	NLayout->addWidget(NSpinBox);
	NLayout->addStretch();
	controlsLayout->addLayout(NLayout);
	connect(NSpinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &MoleculesDynamics::resetAnimation);

	QHBoxLayout* densityLayout = new QHBoxLayout();
	densityLayout->setAlignment(Qt::AlignLeft);
	QLabel* densityLabel = new QLabel("Плотность (ρ)");
	densityLabel->setFixedWidth(180);
	densitySpinBox->setMinimum(0.1);
	densitySpinBox->setMaximum(2);
	densitySpinBox->setSingleStep(0.1);
	densitySpinBox->setValue(0.8);
	densitySpinBox->setFixedWidth(80);
	densityLayout->addWidget(densityLabel);
	densityLayout->addWidget(densitySpinBox);
	densityLayout->addStretch();
	controlsLayout->addLayout(densityLayout);
	connect(densitySpinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &MoleculesDynamics::resetAnimation);

	QHBoxLayout* stepLayout = new QHBoxLayout();
	stepLayout->setAlignment(Qt::AlignLeft);
	QLabel* stepLabel = new QLabel("Шагов моделирования");
	stepLabel->setFixedWidth(180);
	stepSpinBox->setMinimum(100);
	stepSpinBox->setMaximum(100000);
	stepSpinBox->setValue(10000);
	stepSpinBox->setSingleStep(100);
	stepSpinBox->setFixedWidth(80);
	stepLayout->addWidget(stepLabel);
	stepLayout->addWidget(stepSpinBox);
	stepLayout->addStretch();
	controlsLayout->addLayout(stepLayout);
	connect(stepSpinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &MoleculesDynamics::resetAnimation);

	QHBoxLayout* speedLayout = new QHBoxLayout();
	speedLayout->setAlignment(Qt::AlignLeft);
	QLabel* speedLabel = new QLabel("Скорость");
	speedLabel->setFixedWidth(180);
	speedSpinBox->setMinimum(1);
	speedSpinBox->setMaximum(100);
	speedSpinBox->setValue(1);
	speedSpinBox->setFixedWidth(80);
	speedLayout->addWidget(speedLabel);
	speedLayout->addWidget(speedSpinBox);
	speedLayout->addStretch();
	controlsLayout->addLayout(speedLayout);
	connect(speedSpinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &MoleculesDynamics::updateTimerInterval);

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

	for (int i = 0; i < 6; ++i)
	{
		switch (i)
		{
		case 0:
			rungeKutta(method_positions[i], method_velocities[i]);
			break;
		case 1:
			verlet(method_positions[i], method_velocities[i], verlet_prev_positions);
			break;
		case 2:
			velocityVerlet(method_positions[i], method_velocities[i]);
			break;
		case 3:
			leapfrog(method_positions[i], method_velocities[i]);
			break;
		case 4:
			beemanSchofield(method_positions[i], method_velocities[i], method_prev_forces[i]);
			break;
		case 5:
			predictorCorrector(method_positions[i], method_velocities[i], method_prev_forces[i]);
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

void MoleculesDynamics::resetAnimation() {
	currentStep = 0;
	accumulatedMSE.clear();
	accumulatedMSE.resize(6);
	for (double& mse : accumulatedMSE) 
	{
		mse = 0.0;
	}

	int N = NSpinBox->value();
	double L = pow(N / densitySpinBox->value(), 1.0 / 3.0);

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
		double vx = -0.2f + QRandomGenerator::global()->generateDouble() * 0.4f;
		double vy = -0.2f + QRandomGenerator::global()->generateDouble() * 0.4f;
		double vz = -0.2f + QRandomGenerator::global()->generateDouble() * 0.4f;
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
		verlet_prev_positions[j] -= method_velocities[1][j] * 0.0001;
	}

	for (int i = 4; i <= 5; ++i)
	{
		method_prev_forces[i] = computeLJForces(method_positions[i]);
	}

	updateScatters();
	updateMSELabels();
	updateTimerInterval(speedSpinBox->value());
	animationTimer->start();
}

void MoleculesDynamics::resetCamera()
{
	for (int i = 0; i < 6; ++i)
	{
		if (scatters[i])
		{
			Q3DCamera* camera = scatters[i]->scene()->activeCamera();
			camera->setCameraPreset(Q3DCamera::CameraPresetFront);
			camera->setZoomLevel(120);
		}
	}
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

		switch (i)
		{
		case 0:
			rungeKuttaLabel->setText(mseText);
			break;
		case 1:
			verletLabel->setText(mseText);
			break;
		case 2:
			velocityVerletLabel->setText(mseText);
			break;
		case 3:
			leapfrogLabel->setText(mseText);
			break;
		case 4:
			beemanSchofieldLabel->setText(mseText);
			break;
		case 5:
			predictorCorrectorLabel->setText(mseText);
			break;
		}
	}
}

void MoleculesDynamics::updateTimerInterval(double speed)
{
	int interval = static_cast<int>(100 / speed);
	animationTimer->setInterval(interval);
}