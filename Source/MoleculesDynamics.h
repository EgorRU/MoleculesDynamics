#ifndef MOLECULESDYNAMICS_H
#define MOLECULESDYNAMICS_H

#include <random>

#include <QMainWindow>
#include <QGridLayout>
#include <Q3DScatter>
#include <QDoubleSpinBox>
#include <QTimer>
#include <QLabel>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QRandomGenerator>
#include <QMenuBar>
#include <QComboBox>
#include <QPushButton>
#include <QMessageBox>

#include <mathFunc.hpp>

enum class ActionSpinBox 
{
	notGeneration,
	newGeneration,
	notAction,
};

struct SpinBoxConfig
{
	double min;
	double max;
	double step;
	double value;
	QString label;
	double decimals;
	ActionSpinBox method;
};

class MoleculesDynamics : public QMainWindow
{
	Q_OBJECT

public:
	MoleculesDynamics(QWidget* parent = nullptr);

private:
	int currentStep;

	QTimer* animationTimer;
	QGridLayout* gridLayout;
	QLabel* currentTimerLabel;

	QDoubleSpinBox* NSpinBox;
	QDoubleSpinBox* epsilonSpinBox;
	QDoubleSpinBox* sigmaSpinBox;
	QDoubleSpinBox* weightSpinBox;
	QDoubleSpinBox* densitySpinBox;
	QDoubleSpinBox* stepSpinBox;
	QDoubleSpinBox* speedSpinBox;
	QDoubleSpinBox* dtSpinBox;

	QComboBox* MSEComboBox;

	QLabel* rungeKuttaLabel;
	QLabel* verletLabel;
	QLabel* velocityVerletLabel;
	QLabel* leapfrogLabel;
	QLabel* beemanSchofieldLabel;
	QLabel* predictorCorrectorLabel;

	QMap<QDoubleSpinBox*, SpinBoxConfig> spinBoxConfigs;
	QVector<QLabel*> labelsMSE;

	QWidget* centralWidget;
	QWidget* scatterContainers[6];
	Q3DScatter* scatters[6];
	QScatterDataProxy* scatterDataProxies[6];
	QScatter3DSeries* scatterSeries[6];

	QVector<QVector<QVector3D>> method_positions;
	QVector<QVector<QVector3D>> method_velocities;
	QVector<QVector<QVector3D>> method_prev_forces;

	QVector<QVector3D> positions;
	QVector<QVector3D> velocities;
	QVector<QVector3D> verlet_prev_positions;
	QVector<double> accumulatedMSE;

	QVector<QVector3D> initialPositions;
	QVector<QVector3D> initialVelocities;

private:
	void setupUI();
	void start(bool generation);
	void createScatter(int index);
	void updateScatters();
	void positionCorrection(QVector<QVector3D>& pos, QVector<QVector3D>& vel, double L);
	double calculateMSE(const QVector<QVector3D>& pos1, const QVector<QVector3D>& pos2);

private slots:
	void fullRestart();
	void onlyRestart();
	void animateScatters();
	void updateMSELabels();
	void resetCamera();
};

#endif