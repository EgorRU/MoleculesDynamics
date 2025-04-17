#include "MoleculesDynamics.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    QPalette lightPalette;
    lightPalette.setColor(QPalette::Window, QColor(245, 245, 245));
    lightPalette.setColor(QPalette::WindowText, Qt::black);
    lightPalette.setColor(QPalette::Base, Qt::white);
    lightPalette.setColor(QPalette::AlternateBase, QColor(233, 233, 233));
    lightPalette.setColor(QPalette::ToolTipBase, Qt::white);
    lightPalette.setColor(QPalette::ToolTipText, Qt::black);
    lightPalette.setColor(QPalette::Text, Qt::black);
    lightPalette.setColor(QPalette::Button, QColor(240, 240, 240));
    lightPalette.setColor(QPalette::ButtonText, Qt::black);
    lightPalette.setColor(QPalette::BrightText, Qt::red);
    lightPalette.setColor(QPalette::Highlight, QColor(38, 79, 120));
    lightPalette.setColor(QPalette::HighlightedText, Qt::white);

    qApp->setPalette(lightPalette);
    qApp->setStyle("Fusion");

    a.setWindowIcon(QIcon("icons/icon.ico"));
    MoleculesDynamics w;
    w.showMaximized();
    w.setWindowTitle("Движение молекул");
    return a.exec();
}