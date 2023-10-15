#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "equations.h"
#include "atmosphere.h"
#include "airforce.h"
#include "mass.h"
#include "eastwind.h"
#include "block.h"
#include "focus.h"
#include "dycoef.h"
#include <thread>
#include <QtSql>
#include <QDebug>
#include <cmath>
#include <iostream>




MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

// Создание объектов
    block first;
// Исходные данные РН
    double D = 4.1;
    double h = 0.01;

    double M_Rocket = 6545;
    double Peng [2];
    int count;
    double Lmax;
    //double gl_c;
    double gl_cmax;
    double Wind1;
    double Wind2;
    double CILCON;
//
// Глобальные векторы
    QVector<double>lin(0);
    QVector<double>sinn(0);
    QVector<double> jinn(0), jinn2(0);
    QVector<double>CM(0);
    QVector<double>ALI_1(0);
    QVector<double>ALI_2(0);
    QVector<double>res(0);
    QVector<double>TET_1(0);
    QVector<double>TET_2(0);
    QVector<double>center_1(0);
    QVector<double>center_2(0);
    QVector<double>be(0);
    QVector<double>pi(0);
    QVector<double>yu_1(0);
    QVector<double>yu_2(0);
    QVector<double>mass_1(0);
    QVector<double>mass_2(0);
    QVector<double> xn(0);
    QVector<double> v_1(0), v_2(0), v_3(0);
    QVector<double> Long_1(0), Long_2(0), H1(0), H2(0), H3(0);
    QVector<double> angle(0), b1(0), b2(0);
    QVector<double> P1(0), P2(0);
    QVector<double> f1(0);
    QVector<double> Lon(0), Lonre(0);
    QVector<double> w(0), pc2(0);
    QVector<double> cy2(0);
    QVector<double> dyn1(0), dyn2(0);
//
// Кнопка "Трасса полета"
void MainWindow::on_action_triggered()
{


    // Параметры тяги и угла атаки
    QString p1=ui->lineEdit_10->text();
    QString p2=ui->lineEdit_11->text();
    QString al1=ui->lineEdit_5->text();
    QString al2=ui->lineEdit_4->text();

   // Поле локальных переменных
    double Smid = 6.83493;
    first.V = 3750;
    first.tY = 145000;
    first.tX = 450000;
    double alpha;
    double time = 287.7;
    Lmax = 26.681;

    double M_stage[2];
    double T_stage [2];
    double Ott;
    double CX_1, CY_1;
    double CX_2, CY_2;
    count =0;
    double Na=0.1, Mah_1, Mah_2;

    double dN;
    double m_O[2], m_C[2];
    double d_O[2], d_C[2];
    double S_o[2], S_c[2];
    double I_o[2], I_c[2];
    double Imp_t, CF;
    double pitch;
    double bpr;
    double HSP_1, HSP_p_1, VHSP_1;
    double HSP_2, HSP_p_2;
    double MHSP_1;
    double MHSP_2, VHSP_2;
    int e=0;
    double Ott_1, Ott_2 ;
    double H_1;
    int uj=0;
    int ten=0;
    double dep;
    double deo, dec;
    double amax;
    double V1=0, Y1=0, H11=0, X1=0;
    double X_oneC;
    double X_twoC;
    //
    // Создание объектов
    mass M;
    //
    // Определение габаритов ракеты

    //double L = Lmax;
    double L1= M.get_lst1();
    double L2= M.get_lst2();
    //

    // ИД итеративного расчета
        first.m_t = M_Rocket;
        first.anY = 0;
    //
    // Итеративный расчет
    H1.clear();
    H2.clear();
    xn.clear();
    count = 0;

    int i = 0;
    int value = 0;
    d_O[1] = 0;
    //while (second.tY>=0.5 && second.V>=0.5) second.tY>=0.5
    while (first.tY/1000>0.3)
    {
        std::cout << time << " | " << first.V << "  |  " << first.tY/1000 << " | " << first.tX/1000 << std::endl;
        airforce Qus_1 (Mah_1);
        atmosphere Atm_1 (first.tY);
        eastwind W1 (first.tY/1000);
        focus F;
        first.Peng_t = 0;

           // Учет параметров атмосферы

           HSP_p_1 = HSP_1;
           HSP_1 = (Atm_1.get_density()/2)*pow(first.V,2);
           Mah_1 = first.V/Atm_1.get_SV();


           CX_1=Qus_1.getCX();
           CY_1=Qus_1.getCY();

           alpha = 20;
           if (first.tY/1000 >= 20 && first.tY/1000 <= 60)
           {
           CX_1=-0.79;
           CY_1= 1.2;
           alpha = 1;
           }

           equations B_1 (Atm_1.get_density(), Smid, Atm_1.get_AOG(), first.m_t, CX_1, CY_1, first.Peng_t, alpha);

           dN = B_1.fdN(first.tY, first.V, first.anY);
           Ott = first.anY-Na;
           pitch = Ott-first.alpha/57.3;
           first.tX += h * (first.V* cos(first.anY)+X1)/2;
           X1 = first.V* cos(first.anY);
           Wind1 = W1.WSol();

            first.tY += (first.V* sin(first.anY) + H11)/2*h;
            first.V   += (B_1.fdV(first.V, first.anY) + V1)/2*h;
            first.anY += (B_1.fdY(first.tY, first.V, first.anY)+Y1)/2*h;
            V1 = B_1.fdV(first.V, first.anY);
            Y1 = B_1.fdY(first.tY, first.V, first.anY);
            H11 = first.V* sin(first.anY);
            time += h;

                xn.push_back(time);
                yu_1.push_back(HSP_1);
                center_1.push_back(first.gl_c);
                v_1.push_back(first.V);
                sinn.push_back(first.Ssumm);
                jinn.push_back(first.Isumm);
                CM.push_back(first.Ssumm/first.m_t);
                mass_1.push_back(first.m_t);
                Long_1.push_back(first.tX/1000);
                w.push_back(Wind2);
                H1.push_back(first.tY/1000);
                angle.push_back(pitch*57.3);
                b1.push_back(HSP_1*Smid*CX_1);
                TET_1.push_back(first.anY*57.3);
                pi.push_back(pitch*57.3);
                P1.push_back(first.Peng_t);
                f1.push_back(first.focus);
                Lon.push_back(first.Peng_t/(first.m_t*Atm_1.get_AOG()));

        count+=1;
    }
    //db.close();
    QString mk2s = QString::number(first.V);        ui->lineEdit  ->setText(mk2s);
    QString mk1s = QString::number(first.tY/1000);        ui->lineEdit_2->setText(mk1s);




    ui->lineEdit_4->show();
    ui->lineEdit_5->show();
    ui->label_4->show();
    ui->label_5->show();


    ui->action_2->setEnabled(true);


}

void MainWindow::on_action_4_triggered()
{
        close();
}


void MainWindow::on_action_6_triggered()
{
    QVector<double> qu(0), x(0), iner(0);

    for (int i=0;i<=count;i++)

    {
    x.push_back(i*h);
    }
    ui->lineEdit  ->hide();
    ui->lineEdit_2->hide();
    ui->lineEdit_3->hide();
    ui->lineEdit_6->hide();
    ui->lineEdit_7->hide();

    ui->lineEdit_4->hide();
    ui->lineEdit_5->hide();
    ui->label_4->hide();
    ui->label_5->hide();


    ui->label_2->hide();
    ui->label_3->hide();
    ui->label  ->hide();
    ui->label_8->hide();
    ui->label_9->hide();


    ui->label_10->setText("Угол траектории (1), град");
    ui->label_11->setText("Угол траектории (2), град");
}

void MainWindow::on_heightButton_clicked()
{
    ui->widget->clearGraphs();
    ui->widget->addGraph();
    ui->widget->graph(0)->setData(xn,H1);
    ui->widget->xAxis->setLabel("Время, с");
    ui->widget->yAxis->setLabel("Высота, км");
    ui->widget->xAxis->setRange(287.7, count*h+50);
    ui->widget->yAxis->setRange(0, 230);
    ui->widget->replot();

}


void MainWindow::on_distanceButton_clicked()
{
    ui->widget->clearGraphs();
    ui->widget->addGraph();
    ui->widget->graph(0)->setData(xn,Long_1);
    ui->widget->xAxis->setLabel("Время, с");
    ui->widget->yAxis->setLabel("Дальность, км");
    ui->widget->xAxis->setRange(0, count*h);
    ui->widget->yAxis->setRange(0, 1600);
    ui->widget->replot();

    ui->widget->addGraph();
    ui->widget->graph(1)->setData(xn,Long_2);
    ui->widget->graph(1)->setPen(QPen(Qt::PenStyle::DashLine));
    ui->widget->replot();
}


void MainWindow::on_trjButton_clicked()
{
    ui->widget->clearGraphs();
    ui->widget->addGraph();
    ui->widget->graph(0)->setData(Long_1,H1);
    ui->widget->xAxis->setLabel("Дальность, км");
    ui->widget->yAxis->setLabel("Высота, км");
    ui->widget->xAxis->setRange(450, 1500);
    ui->widget->yAxis->setRange(0, 230);
    ui->widget->replot();

}

void MainWindow::on_velButton_clicked()
{
    ui->widget->clearGraphs();
    ui->widget->addGraph();
    ui->widget->graph(0)->setData(xn,v_1);
    ui->widget->xAxis->setLabel("Время, с");
    ui->widget->yAxis->setLabel("Скорость, м/с");
    ui->widget->xAxis->setRange(0, count*h+10);
    ui->widget->yAxis->setRange(0, first.V+10);
    ui->widget->replot();

    ui->widget->addGraph();
    ui->widget->graph(1)->setData(xn,v_2);
    ui->widget->graph(1)->setPen(QPen(Qt::PenStyle::DashLine));
    ui->widget->replot();
}

void MainWindow::on_alphaButton_clicked()
{
    ui->widget->clearGraphs();
    ui->widget->addGraph();
    ui->widget->graph(0)->setData(xn, ALI_1);
    ui->widget->xAxis->setLabel("Время, с");
    ui->widget->yAxis->setLabel("Угол атаки, град");
    ui->widget->xAxis->setRange(0, count*h+20);
    ui->widget->yAxis->setRange(-5, 91);
    ui->widget->replot();

    ui->widget->addGraph();
    ui->widget->graph(1)->setData(xn, ALI_2);
    ui->widget->graph(1)->setPen(QPen(Qt::PenStyle::DashLine));
    ui->widget->replot();
}

void MainWindow::on_angleButton_clicked()
{
    ui->widget->clearGraphs();
    ui->widget->addGraph();
    ui->widget->graph(0)->setData(xn, TET_1);
    ui->widget->graph(0)->setPen(QPen(Qt::red));
    ui->widget->xAxis->setLabel("Время, с");
    ui->widget->yAxis->setLabel("Угол наклона траектории, град");
    ui->widget->xAxis->setRange(0, count*h+10);
    ui->widget->yAxis->setRange(-90, 90);
    ui->widget->replot();

    ui->widget->addGraph();
    ui->widget->graph(1)->setData(xn, TET_2);
    ui->widget->graph(1)->setPen(QPen(Qt::PenStyle::DashLine));
    ui->widget->replot();
}

void MainWindow::on_vhButton_clicked()
{
    ui->widget->clearGraphs();
    ui->widget->addGraph();
    ui->widget->graph(0)->setData(xn, yu_1);
    ui->widget->graph(0)->setPen(QPen(Qt::red));
    ui->widget->xAxis->setLabel("Время, с");
    ui->widget->yAxis->setLabel("Максимальный скоростной напор, Па");
    ui->widget->xAxis->setRange(0, count*h);
    ui->widget->yAxis->setRange(0, *std::max_element(yu_2.begin(),yu_2.end()));
    ui->widget->replot();

    ui->widget->addGraph();
    ui->widget->graph(1)->setData(xn, yu_2);
    ui->widget->graph(1)->setPen(QPen(Qt::PenStyle::DashLine));
    ui->widget->replot();
}

void MainWindow::on_centerButton_clicked()
{
    ui->widget->clearGraphs();
    ui->widget->addGraph();
    ui->widget->graph(0)->setData(xn, center_1);
    ui->widget->graph(0)->setPen(QPen(Qt::green));
    ui->widget->xAxis->setLabel("Время, с");
    ui->widget->yAxis->setLabel("Центр масс РН, м");
    ui->widget->xAxis->setRange(0, count*h);
    ui->widget->yAxis->setRange(0, *std::max_element(center_1.begin(),center_1.end()));
    ui->widget->replot();
}







void MainWindow::on_massButton_clicked()
{
    ui->widget->clearGraphs();
    ui->widget->addGraph();
    ui->widget->graph(0)->setData(xn, mass_1);
    ui->widget->graph(0)->setPen(QPen(Qt::red));
    ui->widget->xAxis->setLabel("Время, с");
    ui->widget->yAxis->setLabel("Масса, кг");
    ui->widget->xAxis->setRange(0, count*h);
    ui->widget->yAxis->setRange(0, *std::max_element(mass_1.begin(),mass_1.end()));
    ui->widget->replot();

    ui->widget->addGraph();
    ui->widget->graph(1)->setData(xn, mass_2);
    ui->widget->graph(1)->setPen(QPen(Qt::PenStyle::DashLine));
    ui->widget->replot();
}


void MainWindow::on_pushButton_11_clicked()
{
    ui->widget->clearGraphs();
    ui->widget->addGraph();
    ui->widget->graph(0)->setData(xn, P1);
    ui->widget->graph(0)->setPen(QPen(Qt::red));
    ui->widget->xAxis->setLabel("Время, с");
    ui->widget->yAxis->setLabel("Тяга ДУ, кН");
    ui->widget->xAxis->setRange(0, count*h);
    ui->widget->yAxis->setRange(0, *std::max_element(P1.begin(),P1.end()));
    ui->widget->replot();

    ui->widget->addGraph();
    ui->widget->graph(1)->setData(xn, P2);
    ui->widget->graph(1)->setPen(QPen(Qt::PenStyle::DashLine));
    ui->widget->replot();
}


void MainWindow::on_pushButton_14_clicked()
{
    ui->widget->clearGraphs();
    ui->widget->addGraph();
    ui->widget->graph(0)->setData(xn, Lon);
    ui->widget->graph(0)->setPen(QPen(Qt::red));
    ui->widget->xAxis->setLabel("Время, с");
    ui->widget->yAxis->setLabel("Продольная перегрузка");
    ui->widget->xAxis->setRange(0, count*h);
    ui->widget->yAxis->setRange(0, 10);
    ui->widget->replot();

    ui->widget->addGraph();
    ui->widget->graph(1)->setData(xn, Lonre);
    ui->widget->graph(1)->setPen(QPen(Qt::PenStyle::DashLine));
    ui->widget->replot();
}


void MainWindow::on_fokusButton_clicked()
{
    ui->widget->clearGraphs();
    ui->widget->addGraph();
    ui->widget->graph(0)->setData(xn, f1);
    ui->widget->graph(0)->setPen(QPen(Qt::red));
    ui->widget->xAxis->setLabel("Время, с");
    ui->widget->yAxis->setLabel("Фокус, м");
    ui->widget->xAxis->setRange(0, count*h);
    ui->widget->yAxis->setRange(0, 5);
    ui->widget->replot();
}


void MainWindow::on_action_5_triggered()
{
    DC d;
    d.data_writing(xn, v_2, H2, w, mass_2, P2, pc2, yu_2, ALI_2, cy2, dyn1, dyn2, jinn2);

    QString k14=ui->lineEdit_14->text(); double P = k14.toDouble();
    QString k15=ui->lineEdit_15->text(); double I = k15.toDouble();
    QString k16=ui->lineEdit_16->text(); double D = k16.toDouble();
    d.get_T().clear();
    d.get_V().clear();
    d.data_calculating(P, I, D);

    ui->widget->clearGraphs();
    ui->widget->addGraph();
    ui->widget->graph(0)->setData(d.get_T(), d.get_V());
    ui->widget->graph(0)->setPen(QPen(Qt::green));
    ui->widget->xAxis->setLabel("Время, с");
    ui->widget->yAxis->setLabel("Отклонение по тангажу, м");
    ui->widget->xAxis->setRange(405, 440);
    ui->widget->yAxis->setRange(-20, 20);
    ui->widget->replot();
        //    query_2.exec("SELECT time, con_thrust, mass, coef_cy, thrust, dy_pressure, X1, inerz, velocity, X2, wind FROM ballistic_data");
}


void MainWindow::on_action_7_triggered()
{

}

