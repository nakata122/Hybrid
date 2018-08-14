#include "mainwindow.h"
#include "main.h"

/**
 * Constructor that connects signals & slots.
 */
MainWindow::MainWindow(QWidget *parent) :
	QMainWindow(parent),
	ui(new Ui::MainWindow) {
	ui->setupUi(this);
    connect(ui->image, SIGNAL(saveImage()),
            this, SLOT(saveImage()));

	connect(this, SIGNAL(smallpaint_fixedSignal(smallpaint_fixed::Vec**, int)),
			this, SLOT(smallpaint_fixedSlot(smallpaint_fixed::Vec**, int)));
	connect(this, SIGNAL(smallpaint_painterlySignal(smallpaint_painterly::Vec**, int)),
			this, SLOT(smallpaint_painterlySlot(smallpaint_painterly::Vec**, int)));
	connect(this, SIGNAL(smallpaint_bvhSignal(smallpaint_bvh::Vec**, int)),
			this, SLOT(smallpaint_bvhSlot(smallpaint_bvh::Vec**, int)));
	connect(this, SIGNAL(smallpaint_pssmltSignal(smallpaint_pssmlt::Vec**, int)),
			this, SLOT(smallpaint_pssmltSlot(smallpaint_pssmlt::Vec**, int)));
	connect(this, SIGNAL(smallpaint_smallmediaSignal(smallpaint_smallmedia::Vec**, int)),
			this, SLOT(smallpaint_smallmediaSlot(smallpaint_smallmedia::Vec**, int)));
	connect(this, SIGNAL(smallpaint_ppmSignal(smallpaint_ppm::Vec**, int)),
			this, SLOT(smallpaint_ppmSlot(smallpaint_ppm::Vec**, int)));
    connect(this, SIGNAL(smallpaint_vrlSignal(smallpaint_vrl::Vec**, int)),
            this, SLOT(smallpaint_vrlSlot(smallpaint_vrl::Vec**, int)));
}

/**
 * Destructor.
 */
MainWindow::~MainWindow() {
	delete ui;
}

/**
 * Displays QMessageBox, closes all windows and cancels rendering for shutdown.
 */
void MainWindow::closeEvent(QCloseEvent *event) {
	QMessageBox::StandardButton resBtn = QMessageBox::question(this, "", tr("Are you sure you want to quit?\nUnsaved renderings will be discarded."), QMessageBox::No | QMessageBox::Yes, QMessageBox::Yes);
	if (resBtn != QMessageBox::Yes) {
		event->ignore();
	} else {
		if (isRendering()) cancelRendering();
		ui->renderInfo->setText("Exit");
		QApplication::closeAllWindows();
		event->accept();
	}
}

/**
 * Centers rendered image in middle of its area.
 */
void MainWindow::resizeEvent(QResizeEvent *event) {
    if (lastImage == "" && image.isNull()) return;
    ui->image->resizeImage(image);
}

void MainWindow::saveImage(){
    if (isRendering()) return;
    helper.saveImage();
}

/**
 * Enter key toggles render button.
 */
void MainWindow::keyPressEvent(QKeyEvent *event) {
	if (event->key() == Qt::Key_Enter || event->key() == Qt::Key_Return) {
		on_renderButton_clicked();
	}
}

/**
 * Handles change of ComboBox.
 */
void MainWindow::on_renderModesComboBox_currentTextChanged(const QString &arg1) {
	helper.changedComboBox(helper.getRenderer(arg1));
	resizeImage();
}

/**
 * Saves new image and displayes it.
 */
void MainWindow::drawImage(QImage img, std::string name, int currentSpp, int goalSpp) {
    image = img.copy();
    if (name == "logo"){
        ui->image->setImage(image);
        resizeImage();
    } else {
        if (!testing) {
            lastImage = name;
            helper.updateInfo(currentSpp, goalSpp);
            ui->image->setImage(image);
            resizeImage();
        } else {
            helper.updateInfo(currentSpp, goalSpp);
            testWindow->setSecondImage(image);
            if (currentSpp == goalSpp) testWindow->setFinished(image);
        }
    }
}

/**
 * Starts or cancels rendering.
 */
void MainWindow::on_renderButton_clicked() {
	if (ui->renderButton->text() == "Render") {
		QString selectedRenderer = helper.getRenderer(ui->renderModesComboBox->currentText());
		if (selectedRenderer == "") return;
		if (selectedRenderer == "tests" || testing) {
			testRender(true);
			return;
		}
		int spp = 0;
		int size = 0;
		float refr = 0;

        smallpaint::RenderInfo info;
        helper.getInput(selectedRenderer, info);
        helper.initializeRenderInformation(selectedRenderer, info);

        drawImage(QImage(":/test_images/logo_processing.png"), "logo", 0, spp);

		start = clock();
        smallpaint::sendToRender(info, "smallpaint_" + selectedRenderer.toStdString());
	} else {
		ui->renderButton->setText("Render");
		if (!isRendering()) {
			ui->renderInfo->setText("Render ended");
			if (testing) {
				on_renderButton_clicked();
				return;
			}
		} else {
			cancelRendering();
			setRendering(false);
			ui->renderInfo->setText("Render canceled");

			if (testing) {
				testing = false;
				testWindow->makeCloseable();
				tests.clear();
			}
		}
		ui->timeInfo->clear();
		ui->progressBar->setValue(0);
		helper.saveImage();
	}
}

/**
 * Resizes the image.
 */
void MainWindow::resizeImage() {
	QResizeEvent* resizeEvent = new QResizeEvent(this->window()->size(), this->window()->size());
	MainWindow::resizeEvent(resizeEvent);
	resizeEvent->~QResizeEvent();
}

/**
 * Handles test renderings.
 */
void MainWindow::testRender(bool renderButtonActivated) {
	if (renderButtonActivated) {
		QString testRenderer;
		int spp = 0;
		int size = 0;
		float refr = 0;

		if (!testing && tests.size() == 0) {
			if (ui->fixedTest->isChecked()) tests.push_back(0);
			if (ui->painterlyTest->isChecked()) tests.push_back(1);
			if (ui->bvhTest->isChecked()) tests.push_back(2);
			if (ui->pssmltTest->isChecked()) tests.push_back(3);
			if (ui->smallmediaTest->isChecked()) tests.push_back(4);
            if (ui->ppmTest->isChecked()) tests.push_back(5);
            if (ui->vrlTest->isChecked()) tests.push_back(6);
        } else if (testing && tests.size() == 0) {
			ui->renderButton->setText("Cancel");
			testing = false;
			tests.clear();
			on_renderButton_clicked();
			return;
		}

        smallpaint::RenderInfo info;
        info.spp = 10000;
        info.size = 512;
        info.refr = 1.5;
		testing = true;

		if (tests.at(0) == 0) {
			testRenderer = "smallpaint_fixed";
		} else if (tests.at(0) == 1) {
			testRenderer = "smallpaint_painterly";
		} else if (tests.at(0) == 2) {
			testRenderer = "smallpaint_bvh";
		} else if (tests.at(0) == 3) {
			testRenderer = "smallpaint_pssmlt";
		} else if (tests.at(0) == 4) {
			testRenderer = "smallpaint_smallmedia";
		} else if (tests.at(0) == 5) {
			testRenderer = "smallpaint_ppm";
            info.spp *= 100000;
        } else if (tests.at(0) == 6) {
            testRenderer = "smallpaint_vrl";
        }
		tests.erase(tests.begin());

        helper.initializeRenderInformation(testRenderer, info);

		testWindow = new TestWindow;
		testWindow->setWindowTitle(testRenderer);
        testWindow->setFirstImage(QImage(":/test_images/" + testRenderer + ".ppm"));
		testWindow->setWindowFlags(Qt::WindowTitleHint);
		testWindow->show();

		start = clock();
        smallpaint::sendToRender(info, testRenderer.toStdString());
	}
}

/**
 * Slots that forward to the processImage functions.
 */

void MainWindow::smallpaint_fixedSlot(smallpaint_fixed::Vec **pix, int spp) {
	smallpaint_fixed::processImage(pix, spp);
}

void MainWindow::smallpaint_painterlySlot(smallpaint_painterly::Vec **pix, int spp) {
	smallpaint_painterly::processImage(pix, spp);
}

void MainWindow::smallpaint_bvhSlot(smallpaint_bvh::Vec **pix, int spp) {
	smallpaint_bvh::processImage(pix, spp);
}

void MainWindow::smallpaint_pssmltSlot(smallpaint_pssmlt::Vec **pix, int spp) {
	smallpaint_pssmlt::processImage(pix, spp);
}

void MainWindow::smallpaint_smallmediaSlot(smallpaint_smallmedia::Vec **pix, int spp) {
	smallpaint_smallmedia::processImage(pix, spp);
}

void MainWindow::smallpaint_ppmSlot(smallpaint_ppm::Vec **pix, int spp) {
	smallpaint_ppm::processImage(pix, spp);
}

void MainWindow::smallpaint_vrlSlot(smallpaint_vrl::Vec **pix, int spp) {
    smallpaint_vrl::processImage(pix, spp);
}

void MainWindow::on_vrlScene_currentIndexChanged(int index) {
    ui->vrlSize->setText(QString::number(512));
    ui->vrlSamples->setText(QString::number(100));
    ui->vrlVRLPS->setText(QString::number(10));
    ui->vrlBounces->setText(QString::number(1));
    ui->vrlRefraction->setText(QString::number(1.5));

    if(index == 0){
        ui->vrlSigmaA->setText(QString::number(.1));
        ui->vrlSigmaS->setText(QString::number(.01));
        ui->vrlG->setText(QString::number(0));
        ui->vrlMediumRadiance->setChecked(true);
        ui->vrlSurfaceRadiance->setChecked(true);
    }else if(index == 1 || index == 2){
        ui->vrlSigmaA->setText(QString::number(.1));
        ui->vrlSigmaS->setText(QString::number(.25));
        ui->vrlG->setText(QString::number(0));
        ui->vrlMediumRadiance->setChecked(true);
        ui->vrlSurfaceRadiance->setChecked(false);
    }else if(index == 3 || index == 4){
        ui->vrlSigmaA->setText(QString::number(.65));
        ui->vrlSigmaS->setText(QString::number(.5));
        ui->vrlG->setText(QString::number(0));
        ui->vrlMediumRadiance->setChecked(index == 3);
        ui->vrlSurfaceRadiance->setChecked(index == 4);
    }else if(index == 5){
        ui->vrlSigmaA->setText(QString::number(.5));
        ui->vrlSigmaS->setText(QString::number(.1));
        ui->vrlG->setText(QString::number(0));
        ui->vrlMediumRadiance->setChecked(true);
        ui->vrlSurfaceRadiance->setChecked(true);
    }
}
