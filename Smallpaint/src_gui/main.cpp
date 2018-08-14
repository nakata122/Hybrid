// smallpaint by karoly zsolnai - zsolnai@cg.tuwien.ac.at
//
// render, modify, create new scenes, tinker around, and most of all:
// have fun!
//
// This program is used as an educational learning tool on the Rendering
// course at TU Wien. Course webpage:
// http://cg.tuwien.ac.at/courses/Rendering/
//
// This UI serves as an educational learning tool to see the progressions
// and experience the different renderes with different input.
// Created by Michael Oppitz - e1227129@student.tuwien.ac.at

#include <QApplication>
#include <QStyleFactory>
#include <mainwindow.h>

void setStyleSheet(){
    qApp->setStyle(QStyleFactory::create("Fusion"));

    QPalette darkPalette;
    darkPalette.setColor(QPalette::Window, QColor(53,53,53));
    darkPalette.setColor(QPalette::WindowText, Qt::white);
    darkPalette.setColor(QPalette::Base, QColor(25,25,25));
    darkPalette.setColor(QPalette::AlternateBase, QColor(53,53,53));
    darkPalette.setColor(QPalette::ToolTipBase, Qt::white);
    darkPalette.setColor(QPalette::ToolTipText, Qt::white);
    darkPalette.setColor(QPalette::Text, Qt::white);
    darkPalette.setColor(QPalette::Button, QColor(21, 65, 109));
    darkPalette.setColor(QPalette::ButtonText, Qt::white);
    darkPalette.setColor(QPalette::BrightText, Qt::red);
    darkPalette.setColor(QPalette::Link, QColor(21, 65, 109));
    darkPalette.setColor(QPalette::Shadow, QColor(25, 25, 25));
    darkPalette.setColor(QPalette::Midlight, QColor(53, 53, 53));
    darkPalette.setColor(QPalette::Dark, QColor(25, 25, 25));
    darkPalette.setColor(QPalette::Light, QColor(75, 75, 75));

    darkPalette.setColor(QPalette::Highlight, QColor(21, 65, 109));
    darkPalette.setColor(QPalette::HighlightedText, Qt::black);

    qApp->setPalette(darkPalette);
}

/**
 * Creates new MainWindow and sets instances.
 */
int main(int argc, char *argv[]) {
    QApplication a(argc, argv);
    setStyleSheet();

	MainWindow *w = new MainWindow;
	setInstanceForRenderDistribution(w);
	setInstanceForHelperFunctions(w);
	w->show();
	return a.exec();
}
