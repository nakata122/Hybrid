#include "image.h"
#include <qdebug.h>
Image::Image(QWidget *parent) : QWidget(parent) {
    QHBoxLayout* layout = new QHBoxLayout;

    _label = new QLabel;
    _label->setAlignment(Qt::AlignCenter);
    _label->setPixmap(QPixmap::fromImage(QImage(":/test_images/logo.png")).scaled(512, 512));
    layout->addWidget(_label, Qt::AlignCenter);
    this->setLayout(layout);
}

void Image::setImage(QImage img) {
    _label->setPixmap(QPixmap::fromImage(img));
}

void Image::resizeImage(QImage img) {
    float value = std::min(this->size().width(), this->size().height());
    if (img.size().width() != 0) _label->setPixmap(QPixmap::fromImage(img.scaled(value, value)));
}

void Image::mouseReleaseEvent(QMouseEvent *event) {
    if (event->button() == Qt::RightButton)
        emit this->saveImage();
}
