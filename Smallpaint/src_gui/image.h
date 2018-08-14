#ifndef IMAGE_H
#define IMAGE_H

#include <QWidget>
#include <QImage>
#include <QLabel>
#include <QHBoxLayout>
#include <QMouseEvent>

class Image : public QWidget
{
    Q_OBJECT
public:
    explicit Image(QWidget *parent = nullptr);
    void setImage(QImage img);
    void resizeImage(QImage img);

signals:
    void saveImage();

public slots:

protected:
    void mouseReleaseEvent(QMouseEvent *event) override;

private:
    QLabel *_label;
};

#endif // IMAGE_H
