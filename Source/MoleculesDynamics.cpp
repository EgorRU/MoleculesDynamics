#include "MoleculesDynamics.h"

QString labels[6] =
{
	"Метод Рунге-Кутта 4 порядка",
	"Метод Верле",
	"Скоростной метод Верле",
	"Метод с перескоками",
	"Метод Бимана-Шофилда",
	"Метод предиктор-корректор"
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
	dtSpinBox = new QDoubleSpinBox();

	MSEComboBox = new QComboBox();

	rungeKuttaLabel = new QLabel();
	verletLabel = new QLabel();
	velocityVerletLabel = new QLabel();
	leapfrogLabel = new QLabel();
	beemanSchofieldLabel = new QLabel();
	predictorCorrectorLabel = new QLabel();

	spinBoxConfigs =
	{
		{stepSpinBox, {100, 100000, 100, 10000, "Шагов моделирования", 0, ActionSpinBox::notAction}},
		{NSpinBox, {2, 1000, 1, 10, "Количество молекул", 0, ActionSpinBox::newGeneration}},
		{densitySpinBox, {0.1, 15, 0.1, 5, "Плотность молекул (ρ)", 2, ActionSpinBox::notGeneration}},
		{weightSpinBox, {0.1, 1000, 0.1, 0.1, "Масса", 2, ActionSpinBox::notGeneration}},
		{epsilonSpinBox, {0.1, 100, 0.05, 1, "ε", 2, ActionSpinBox::notGeneration}},
		{sigmaSpinBox, {0.1, 100, 0.05, 1, "σ", 2, ActionSpinBox::notGeneration}},
		{speedSpinBox, {0.05, 10, 0.05, 2, "Δ начальных скоростей", 1, ActionSpinBox::notAction}},
		{dtSpinBox, {0.000001, 0.001, 0.000001, 0.000001, "Шаг интегрирования (dt)", 6, ActionSpinBox::notGeneration}},
	};

	labelsMSE = QVector<QLabel*>
	{
		rungeKuttaLabel,
		verletLabel,
		velocityVerletLabel,
		leapfrogLabel,
		beemanSchofieldLabel,
		predictorCorrectorLabel,
	};

	for (auto label : labelsMSE)
	{
		label->setTextInteractionFlags
		(
			Qt::TextSelectableByMouse | Qt::TextSelectableByKeyboard
		);
	}

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
		spinBox->setDecimals(config.decimals);
		ActionSpinBox method = config.method;

		QHBoxLayout* layout = new QHBoxLayout();
		layout->setAlignment(Qt::AlignLeft);
		QLabel* label = new QLabel(config.label);
		label->setFixedWidth(180);
		layout->addWidget(label);
		layout->addWidget(spinBox);
		layout->addStretch();

		QString generationText;
		switch (method) {
		case ActionSpinBox::newGeneration:
			generationText = "Генерация+перезапуск";
			break;
		case ActionSpinBox::notGeneration:
			generationText = "Перезапуск";
			break;
		case ActionSpinBox::notAction:
			generationText = "Без действия";
			break;
		}
		QLabel* generationLabel = new QLabel(generationText);
		generationLabel->setFixedWidth(150);
		layout->addWidget(generationLabel);
		controlsLayout->addLayout(layout);

		switch (method)
		{
		case ActionSpinBox::newGeneration:
			connect(spinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
				this, &MoleculesDynamics::fullRestart);
			break;

		case ActionSpinBox::notGeneration:
			connect(spinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
				this, &MoleculesDynamics::onlyRestart);
			break;

		case ActionSpinBox::notAction:
			break;
		}
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
	connect(MSEComboBox, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &MoleculesDynamics::onlyRestart);

	controlsLayout->addWidget(rungeKuttaLabel);
	controlsLayout->addWidget(verletLabel);
	controlsLayout->addWidget(velocityVerletLabel);
	controlsLayout->addWidget(leapfrogLabel);
	controlsLayout->addWidget(beemanSchofieldLabel);
	controlsLayout->addWidget(predictorCorrectorLabel);

	QPushButton* regenerateButton = new QPushButton("Сгенерировать новые начальные условия");
	regenerateButton->setFixedWidth(400);
	controlsLayout->addWidget(regenerateButton);
	connect(regenerateButton, &QPushButton::clicked, this, &MoleculesDynamics::fullRestart);

	QPushButton* restartButton = new QPushButton("Начать моделирование заново на старых начальных условиях");
	restartButton->setFixedWidth(400);
	controlsLayout->addWidget(restartButton);
	connect(restartButton, &QPushButton::clicked, this, &MoleculesDynamics::onlyRestart);

	QHBoxLayout* buttonLayout = new QHBoxLayout();
	buttonLayout->setAlignment(Qt::AlignLeft);

	QPushButton* pauseButton = new QPushButton("Пауза/Продолжить");
	pauseButton->setFixedWidth(200);
	buttonLayout->addWidget(pauseButton);
	connect(pauseButton, &QPushButton::clicked, this, [this]()
		{
			if (animationTimer->isActive())
			{
				animationTimer->stop();
			}
			else
			{
				animationTimer->start();
			}
		});

	QPushButton* nextStepButton = new QPushButton("Следующий шаг");
	nextStepButton->setFixedWidth(200);
	buttonLayout->addWidget(nextStepButton);
	connect(nextStepButton, &QPushButton::clicked, this, [this]()
		{
			if (!animationTimer->isActive())
			{
				animateScatters();
			}
		});

	controlsLayout->addLayout(buttonLayout);

	QPushButton* cameraButton = new QPushButton("Сброс камер для графиков");
	cameraButton->setFixedWidth(400);
	controlsLayout->addWidget(cameraButton);
	connect(cameraButton, &QPushButton::clicked, this, &MoleculesDynamics::resetCamera);

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
	fullRestart();
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

		animationTimer->stop();

		QString mseMessage = "Итоговые значения MSE:\n\n";
		int referenceMethod = MSEComboBox->currentIndex();

		for (int i = 0; i < 6; ++i)
		{
			double averageMSE = accumulatedMSE[i] / totalSteps;
			mseMessage += QString("%1: %2\n").arg(labels[i]).arg(averageMSE, 0, 'f', 10);
		}

		QMessageBox msgBox;
		msgBox.setWindowTitle("Моделирование завершено");
		msgBox.setText("Моделирование достигло максимального количества шагов");
		msgBox.setInformativeText(mseMessage);
		msgBox.setStandardButtons(QMessageBox::Ok);

		if (QLabel* label = msgBox.findChild<QLabel*>()) 
		{
			label->setTextInteractionFlags(Qt::TextSelectableByMouse | Qt::TextSelectableByKeyboard);
			label->setCursor(Qt::IBeamCursor);
		}

		msgBox.exec();
		return;
	}

	double epsilon = epsilonSpinBox->value();
	double sigma = sigmaSpinBox->value();
	double weight = weightSpinBox->value();
	double dt = dtSpinBox->value();

	for (int i = 0; i < 6; ++i)
	{
		switch (i)
		{
		case 0:
			rungeKutta(method_positions[i], method_velocities[i], epsilon, sigma, weight, dt);
			break;
		case 1:
			verlet(method_positions[i], method_velocities[i], verlet_prev_positions, epsilon, sigma, weight, dt);
			break;
		case 2:
			velocityVerlet(method_positions[i], method_velocities[i], epsilon, sigma, weight, dt);
			break;
		case 3:
			leapfrog(method_positions[i], method_velocities[i], epsilon, sigma, weight, dt);
			break;
		case 4:
			beemanSchofield(method_positions[i], method_velocities[i], method_prev_forces[i], epsilon, sigma, weight, dt);
			break;
		case 5:
			predictorCorrector(method_positions[i], method_velocities[i], method_prev_forces[i], epsilon, sigma, weight, dt);
			break;
		}
		double L = pow(NSpinBox->value() / densitySpinBox->value(), 1.0 / 3.0);
		positionCorrection(method_positions[i], method_velocities[i], L);
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

void MoleculesDynamics::fullRestart()
{
	bool generation = true;
	start(generation);
}

void MoleculesDynamics::onlyRestart()
{
	bool generation = false;
	start(generation);
}

void MoleculesDynamics::start(bool generation)
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
	double epsilon = epsilonSpinBox->value();
	double sigma = sigmaSpinBox->value();
	double weight = weightSpinBox->value();
	double dt = dtSpinBox->value();

	// новая генерация
	if (generation)
	{
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

		static std::mt19937 gen(std::random_device{}());
		static std::normal_distribution<double> dist(0.0, speed / 3.0);

		for (int i = 0; i < N; ++i)
		{
			velocities[i] = QVector3D(dist(gen), dist(gen), dist(gen));
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

		initialPositions = positions;
		initialVelocities = velocities;
	}

	positions = initialPositions;
	velocities = initialVelocities;

	method_positions.clear();
	method_velocities.clear();
	method_prev_forces.clear();

	for (int i = 0; i < 6; ++i)
	{
		method_positions.push_back(positions);
		method_velocities.push_back(velocities);
		method_prev_forces.push_back(QVector<QVector3D>(N, QVector3D(0, 0, 0)));
	}

	QVector<QVector3D> init_forces = computeLJForces(method_positions[1], epsilon, sigma, weight);
	verlet_prev_positions = method_positions[1];
	for (int j = 0; j < N; ++j)
	{
		QVector3D acceleration = init_forces[j] / weight;
		verlet_prev_positions[j] -= method_velocities[1][j] * dt;
		verlet_prev_positions[j] += 0.5 * acceleration * dt * dt;
	}

	for (int i = 4; i <= 5; ++i)
	{
		method_prev_forces[i] = computeLJForces(method_positions[i], epsilon, sigma, weight);
	}

	updateScatters();
	updateMSELabels();
	animationTimer->start();
}

void MoleculesDynamics::positionCorrection(QVector<QVector3D>& pos, QVector<QVector3D>& vel, double L)
{
	for (int i = 0; i < pos.size(); ++i)
	{
		double vx = vel[i].x();
		double vy = vel[i].y();
		double vz = vel[i].z();

		double x = reflect(pos[i].x(), L, vx);
		double y = reflect(pos[i].y(), L, vy);
		double z = reflect(pos[i].z(), L, vz);

		pos[i].setX(x);
		pos[i].setY(y);
		pos[i].setZ(z);

		vel[i].setX(vx);
		vel[i].setY(vy);
		vel[i].setZ(vz);
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
		double averageMSE = (currentStep > 0) ? accumulatedMSE[i] / currentStep : 0.0;
		QString mseText = QString("MSE %1: %2").arg(labels[i]).arg(averageMSE, 0, 'f', 10);
		labelsMSE[i]->setText(mseText);
	}
}

void MoleculesDynamics::resetCamera()
{
	for (int i = 0; i < 6; ++i)
	{
		auto cam = scatters[i]->scene()->activeCamera();
		cam->setCameraPreset(Q3DCamera::CameraPresetFront);
		cam->setZoomLevel(120);
	}
}