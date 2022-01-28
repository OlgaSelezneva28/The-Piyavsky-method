using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using ZedGraph;


namespace Piyavsky_Method
{

    public partial class Form1 : Form
    {
        int N;
        int M;
        int number; // номер тестовой задачи 

        public Form1()
        {
            InitializeComponent();

            N = 2;
            double E = double.Parse(textBox4.Text);
            //n (0.5) ^2m <= E
            int mm = 0; // Желательное M
            while ((N * Math.Pow(0.5, 2.0 * mm) >= E) && (N * mm < 52))
            {
                mm++;
            }

            textBox8.Text = Convert.ToString(mm);
        }


        //Функции
        double F(double x1, double x2)
        {
            if (number == 1)
                return Ackley(x1, x2);
            if (number == 2)
                return Levi13(x1, x2);
            if (number == 3)
                return Bukina6(x1, x2);

            textBox1.Text += " Нет такого примера ";
            return -1.0;
        }
        //1
        double Ackley(double x, double y)
        {
            return (-20.0 * Math.Exp(-0.2 * Math.Sqrt(0.5 * (Math.Pow(x, 2) + Math.Pow(y, 2)))) - Math.Exp( 0.5 * (Math.Cos(2.0 * Math.PI * x) + Math.Cos(2 * Math.PI * y))) + Math.Exp(1) + 20.0 );
        }

        //2
        double Levi13(double x, double y)
        {
            return (Math.Pow(Math.Sin(3 * Math.PI * x), 2) + Math.Pow(x - 1.0, 2) * (1 + Math.Pow(Math.Sin(3 * Math.PI * y), 2)) + Math.Pow(y - 1, 2) * (1 + Math.Pow(Math.Sin(3 * Math.PI * y), 2)));
        }


        //3
        double Bukina6(double x, double y)
        {
            return (100.0 * Math.Sqrt(Math.Abs(y - 0.01 * Math.Pow(x, 2))) + 0.01 * Math.Abs( x + 10 ));
        }

        //F(l(x))
        double f(double x)
        {
            /*
            //Для каждой точки ищем центр
            double[] y = new double[2];
            y = Piyavsky_Method.PeanoMap.CalculateImage(x, N, M);

            // Значение исходной функции в данном узле 
            return F(y[0], y[1]);    
             * */
            return Math.Sin(2 * x) + 3 * Math.Sin(3 * x);
        }

         
        //Оценка константы Липшица
        double LK(double r, List<double> X)
        {
            double MuMax = Math.Abs(f(X[1]) - f(X[0])) / (X[1] - X[0]);

            for (int i = 1; i < X.Count - 1; i++)
            {
                double Mu = Math.Abs(f(X[i + 1]) - f(X[i])) / (X[i + 1] - X[i]);
                if (Mu > MuMax)
                    MuMax = Mu;
            }

            //Оценка
            double mp;
            if (MuMax == 0)
                mp = 1;
            else
                mp = r * MuMax;

            return mp;
        }

        
    
        //Запуск
        private void button1_Click(object sender, EventArgs e)
        {
            /*
            double left_border1 = double.Parse(textBox6.Text);
            double left_border2 = double.Parse(textBox7.Text);
            double right_border1 = double.Parse(textBox2.Text);
            double right_border2 = double.Parse(textBox3.Text);
            int razbienia = int.Parse(textBox10.Text);
            number = int.Parse(textBox9.Text); 

            //Очистка строк и столбцов таблицы
            dataGridView1.Rows.Clear();
            dataGridView1.Columns.Clear();

             N = 2;
             M = int.Parse(textBox8.Text);


            // Метод Пиявского
            double E = double.Parse(textBox4.Text);
            double rr = double.Parse(textBox1.Text); // Параметр надежности
            double max_step = int.Parse(textBox5.Text);  // Количество итераций

            List<double> f2_list = new List<double>(); //Точки испытаний

            List<double> X = new List<double>();

            double left_border = 0;
            double right_border = 1;

            f2_list.Add(left_border);
            f2_list.Add(right_border);

            X.Add(left_border);
            X.Add(right_border);

            //  Для выходных данных
            double x_min = left_border;  //  содержит текущий минимум
            int i_min = 0;  //  Содержит номер интервала, на котором находится текущий минимум
            double x_max = left_border;

            
            for (double i = left_border; i <= right_border; i += (right_border - left_border) / razbienia)
            {
                if (f(x_max) < f(i))
                {
                    x_max = i;
                }
            }


            //Проверяем граничные значения
            if (f(x_min) > f(right_border))
                x_min = right_border;

            if (f(x_max) < f(right_border))
                x_max = right_border;

            int t; //  номер интервала максимальной характеристики
            double eps;
            int p = 0;
            double L = rr;

            do
            {
                //Оценка константы Липшица
                L = LK(rr, X);
                //

                List<double> R = new List<double>();

                //  Вычислили характеристики
                for (int i = 0; i < X.Count - 1; i++)
                {
                    R.Add(L * (X[i + 1] - X[i]) / 2.0 - (f(X[i + 1]) + f(X[i])) / 2.0);
                }

                double Rmax = R[0];
                t = 0;
                for (int i = 0; i < R.Count - 1; i++)
                {
                    if (R[i + 1] > Rmax)
                    {
                        t = i + 1;
                        Rmax = R[i + 1];
                    }
                }

                //  Вычисление новой точки
                double x_new = ((X[t + 1] + X[t]) / 2.0) - ((f(X[t + 1]) - f(X[t])) / (2.0 * L));
                if (f(x_min) > f(x_new)) // сравнение Qk* and Qk+1
                {
                    x_min = x_new;
                    i_min = t;
                }

                // Qk* - Rt
                eps = X[i_min + 1] - X[i_min];

                f2_list.Add(x_new);

                X.Add(x_new);
                X.Sort();
                p++;

            } while ((p < max_step) && (eps >= E));

            double U = f(x_min) - Math.Pow(0.5, (M + 1)) * L * Math.Sqrt(N); // нижняя оценка F(y) - исходная задача 
            //Найдем координаты центра гиперкуба для полученной точки 
            double[] yi = new double[2];
            yi = Piyavsky_Method.PeanoMap.CalculateImage(x_min, N, M);
            //Найдем исходные точки 
            double x1 = yi[0] * (right_border1 - left_border1) + (right_border1 + left_border1) / 2.0;
            double x2 = yi[1] * (right_border2 - left_border2) + (right_border2 + left_border2) / 2.0;


            linkLabel1.Text = Convert.ToString(eps);
            linkLabel2.Text = Convert.ToString(p);

            linkLabel3.Text = Convert.ToString(x1);
            linkLabel4.Text = Convert.ToString(x2);
            linkLabel5.Text = Convert.ToString(U);


            //Точки проведения испытаний метода Пиявского
            dataGridView1.RowCount = f2_list.Count + 1;
            dataGridView1.ColumnCount = 4;

            dataGridView1.Rows[0].Cells[0].Value = "№";
            dataGridView1.Rows[0].Cells[1].Value = "0 << xi << 1";
            dataGridView1.Rows[0].Cells[2].Value = "F(l(xi))";

            for (int i = 0; i < f2_list.Count; i++)
            {
                dataGridView1.Rows[i + 1].Cells[0].Value = i;
                dataGridView1.Rows[i + 1].Cells[1].Value = f2_list[i];
                dataGridView1.Rows[i + 1].Cells[2].Value = f(f2_list[i]);
            }

            */

            GraphPane panel = zedGraphControl1.GraphPane;
	        panel.CurveList.Clear();
	        PointPairList f1_list = new PointPairList();
	        PointPairList f2_list = new PointPairList();
	        PointPairList min_point = new PointPairList();

            double left_border = double.Parse(textBox6.Text);
            double right_border = double.Parse(textBox2.Text);

            // Метод Пиявского
            double E = double.Parse(textBox4.Text);
            double rr = double.Parse(textBox1.Text); // Параметр надежности
            double max_step = int.Parse(textBox5.Text);  // Количество итераций

           

            List<double> X = new List<double>();

            f2_list.Add(left_border, 0);
            f2_list.Add(right_border, 0);

            X.Add(left_border);
            X.Add(right_border);

            //  Для выходных данных
            double x_min = left_border;  //  содержит текущий минимум
            int i_min = 0;  //  Содержит номер интервала, на котором находится текущий минимум
            double x_max = left_border;


            for (double i = left_border; i < right_border; i += (right_border - left_border) / 100)
            {
                if (f(x_max) < f(i))
                {
                    x_max = i;
                }
                f1_list.Add(i, f(i));
            }
            f1_list.Add(right_border, f(right_border));

            //Проверяем граничные значения
            if (f(x_min) > f(right_border))
                x_min = right_border;

            if (f(x_max) < f(right_border))
                x_max = right_border;

            int t; //  номер интервала максимальной характеристики
            double eps;
            int p = 0;
            double L = rr;

            do
            {
                //Оценка константы Липшица
                L = LK(rr, X);
                //

                List<double> R = new List<double>();

                //  Вычислили характеристики
                for (int i = 0; i < X.Count - 1; i++)
                {
                    R.Add(L * (X[i + 1] - X[i]) / 2.0 - (f(X[i + 1]) + f(X[i])) / 2.0);
                }

                double Rmax = R[0];
                t = 0;
                for (int i = 0; i < R.Count - 1; i++)
                {
                    if (R[i + 1] > Rmax)
                    {
                        t = i + 1;
                        Rmax = R[i + 1];
                    }
                }

                //  Вычисление новой точки
                double x_new = ((X[t + 1] + X[t]) / 2.0) - ((f(X[t + 1]) - f(X[t])) / (2.0 * L));
                if (f(x_min) > f(x_new)) // сравнение Qk* and Qk+1
                {
                    x_min = x_new;
                    i_min = t;
                }

                // Qk* - Rt
                eps = X[i_min + 1] - X[i_min];

                f2_list.Add(x_new, 0);

                X.Add(x_new);
                X.Sort();
                p++;

            } while ((p < max_step) && (eps >= E));
            min_point.Add(x_min, f(x_min));



            // Устанавливаем интересующий нас интервал по оси X
            panel.XAxis.Scale.Min = left_border - 1;
            panel.XAxis.Scale.Max = right_border + 1;

            // Устанавливаем интересующий нас интервал по оси Y
            panel.YAxis.Scale.Min = f(left_border) - 1;
            panel.YAxis.Scale.Max = f(right_border) + 1;

            LineItem Curve1 = panel.AddCurve("График функции", f1_list, Color.Red, SymbolType.None);
            LineItem Curve2 = panel.AddCurve("Точки испытаний", f2_list, Color.Blue, SymbolType.Circle);
            LineItem Curve3 = panel.AddCurve("Глобальный минимум", min_point, Color.Green, SymbolType.Star);

            zedGraphControl1.AxisChange();
            // Обновляем график
            zedGraphControl1.Invalidate();
        }


        private void zedGraphControl1_Load(object sender, EventArgs e)
        {

        }

        private void button2_Click(object sender, EventArgs e)
        {

            //Очистка строк и столбцов таблицы
            dataGridView1.Rows.Clear();
            dataGridView1.Columns.Clear();


        }

        private void label2_Click(object sender, EventArgs e)
        {

        }

        //Посчитать M
        private void button3_Click(object sender, EventArgs e)
        {
            N = 2;
            double E = double.Parse(textBox4.Text);
            //n (0.5) ^2m <= E
            int mm = 0; // Желательное M
            while ((N * Math.Pow(0.5, 2.0 * mm) >= E) && (N * mm < 52))
            {
                mm++;
            }

            textBox8.Text = Convert.ToString(mm);
        }

        private void Form1_Load(object sender, EventArgs e)
        {

        }

        //Параметры для тестовых функций
        private void button4_Click(object sender, EventArgs e)
        {
            number = int.Parse(textBox9.Text);

            if (number == 1)//Экли
            {
                textBox6.Text = Convert.ToString("-5");
                textBox2.Text = Convert.ToString("5");
                textBox7.Text = Convert.ToString("-5");
                textBox3.Text = Convert.ToString("5");
            }

            if (number == 2)//Леви
            {
                textBox6.Text = Convert.ToString("-10");
                textBox2.Text = Convert.ToString("10");
                textBox7.Text = Convert.ToString("-10");
                textBox3.Text = Convert.ToString("10");
            }

            if (number == 3)//Букина 6
            {
                textBox6.Text = Convert.ToString("-15");
                textBox2.Text = Convert.ToString("-5");
                textBox7.Text = Convert.ToString("-3");
                textBox3.Text = Convert.ToString("3");
            }
        }

    }
}
