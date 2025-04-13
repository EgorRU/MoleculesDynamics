#ifndef MOLECULESDYNAMICS_H
#define MOLECULESDYNAMICS_H

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

#include <mathFunc.hpp>

class MoleculesDynamics : public QMainWindow
{
	Q_OBJECT

public:
	MoleculesDynamics(QWidget* parent = nullptr);

private:
	int currentStep;
	int animationSpeed;

	QTimer* animationTimer;
	QGridLayout* gridLayout;
	QWidget* centralWidget;
	QWidget* scatterContainers[6];
	Q3DScatter* scatters[6];
	QScatterDataProxy* scatterDataProxies[6];
	QScatter3DSeries* scatterSeries[6];
	QComboBox* MSEComboBox;

	QDoubleSpinBox* NSpinBox;
	QDoubleSpinBox* densitySpinBox;
	QDoubleSpinBox* stepSpinBox;
	QDoubleSpinBox* speedSpinBox;

	QLabel* currentTimerLabel;
	QLabel* speedLabel;
	QLabel* rungeKuttaLabel;
	QLabel* verletLabel;
	QLabel* velocityVerletLabel;
	QLabel* leapfrogLabel;
	QLabel* beemanSchofieldLabel;
	QLabel* predictorCorrectorLabel;

	QVector<QVector<QVector3D>> method_positions;
	QVector<QVector<QVector3D>> method_velocities;
	QVector<QVector<QVector3D>> method_prev_forces;

	QVector<QVector3D> positions;
	QVector<QVector3D> velocities;
	QVector<QVector3D> verlet_prev_positions;
	QVector<double> accumulatedMSE;

private:
	void setupUI();
	void createScatter(int index);
	void updateScatters();
	void resetAnimation();
	void applyBoundaryConditions(QVector<QVector3D>& pos, QVector<QVector3D>& vel, double L);
	double calculateMSE(const QVector<QVector3D>& pos1, const QVector<QVector3D>& pos2);

private slots:
	void resetCamera();
	void animateScatters();
	void updateMSELabels();
	void updateTimerInterval(double speed);
};

#endif