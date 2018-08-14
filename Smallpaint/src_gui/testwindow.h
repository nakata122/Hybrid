#ifndef TESTWINDOW_H
#define TESTWINDOW_H

#include <QDialog>
#include <QCloseEvent>
#include <QMessageBox>
#include <QFileDialog>

namespace Ui {
class TestWindow;
}

class TestWindow : public QDialog {
	Q_OBJECT

public:
	explicit TestWindow(QWidget *parent = 0);
	~TestWindow();
	void setFirstImage(QImage img);
	void setSecondImage(QImage img);
	void makeCloseable();
	double analyzeImages();
	void setFinished(QImage img);

protected:
	void mouseReleaseEvent(QMouseEvent *event);

private:
	Ui::TestWindow *ui;
	bool finished = false;
	QImage comparisonImage;
	QImage testImage;
};

#endif // TESTWINDOW_H
