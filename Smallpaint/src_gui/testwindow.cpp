#include "testwindow.h"
#include "ui_testwindow.h"

TestWindow::TestWindow(QWidget *parent) :
	QDialog(parent),
	ui(new Ui::TestWindow) {
	ui->setupUi(this);
}

TestWindow::~TestWindow() {
	delete ui;
}

void TestWindow::setFirstImage(QImage img) {
	comparisonImage = img;
	//ui->label->setText("Rendering...");
    ui->img1->setPixmap(QPixmap::fromImage(img).scaled(512, 512));
}

void TestWindow::setSecondImage(QImage img) {
	ui->img2->setPixmap(QPixmap::fromImage(img));
}

void TestWindow::makeCloseable() {
	QPoint pos = this->pos();
	this->setWindowFlags(Qt::WindowCloseButtonHint);
	this->show();
	this->move(pos);
}

double TestWindow::analyzeImages() {
	int error = 0;
	for (int row = 0; row < testImage.height(); row++) {
		for (int col = 0; col < testImage.width(); col++) {
			QColor testColor(testImage.pixel(row, col));
			int testR = testColor.red();
			int testG = testColor.green();
			int testB = testColor.blue();

			QColor comparisonColor(comparisonImage.pixel(row, col));
			int comparisonR = comparisonColor.red();
			int comparisonG = comparisonColor.green();
			int comparisonB = comparisonColor.blue();

			int threshold = 1;

			if (!(comparisonR - threshold <= testR && testR <= comparisonR + threshold) ||
				!(comparisonG - threshold <= testG && testG <= comparisonG + threshold) ||
				!(comparisonB - threshold <= testB && testB <= comparisonB + threshold)) error++;
		}
	}
	return (double)error / (double)(testImage.height()*testImage.width());
}

void TestWindow::setFinished(QImage img) {
	finished = true;
	testImage = img;
	makeCloseable();
	//double errorRate = analyzeImages() * 100;
	//ui->label->setText("Error Rate: " + QString::number(errorRate, 'f', 2) + "%");
}

void TestWindow::mouseReleaseEvent(QMouseEvent *event) {
	if (!finished) return;
	if (event->button() == Qt::RightButton) {
		if (event->pos().x() >= comparisonImage.width() && event->pos().y() <= comparisonImage.height()) {
			QString filename = QFileDialog::getSaveFileName(this, tr("Save Image"), this->windowTitle(), tr("Image Files (*.ppm)"));
			testImage.save(filename, "PPM", -1);
		}
	}
}
