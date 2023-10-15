#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include <QMainWindow>
#include <qcustomplot.h>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:


    void on_action_triggered();

    void on_action_4_triggered();

    void on_action_6_triggered();

    void on_heightButton_clicked();

    void on_distanceButton_clicked();

    void on_trjButton_clicked();

    void on_velButton_clicked();

    void on_alphaButton_clicked();

    void on_angleButton_clicked();

    void on_vhButton_clicked();

    void on_centerButton_clicked();

    void on_massButton_clicked();

    void on_pushButton_11_clicked();

    void on_pushButton_14_clicked();

    void on_fokusButton_clicked();

    void on_action_5_triggered();

    void on_action_7_triggered();

private:
    Ui::MainWindow *ui;
    QCustomPlot *wGraphic;
    QCPCurve *verticalLine;
    QCPItemTracer *tracer;
};
#endif // MAINWINDOW_H
