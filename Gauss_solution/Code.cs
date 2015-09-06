//Добрий день! Відповідно до завдання, матриця коефіцієнтів за исана в коді. 
//Код спопіювати і вставити в новий проект Visual Stusio. 
// Натистнути Ctrl+F5. Результат обчислння виводиться в консолі.

//Завчасно дякую за Вашу увагу!



using System;
using System.Linq;

namespace Gauss
{
    class Program
    {
        static void Main()
        {
            double[,] matrix = {
	                {-3, 4, 1, 4},
	                {0, 1, 3, 2},
	                {4, 0, -2, -3},
	                {1000, 3, 1, -5} 
	                               }; // Коефіцієнти системи рівнянь
            double[] b = { -1, -1, 4, -2 }; // Права частина системи рівнень.


            try
            {
                LinearSystem ls = new LinearSystem(matrix, b, 0.0001);
                Console.WriteLine("Вектор розв'язку Х: ");
                Console.WriteLine(String.Join("\n", ls.XVector));
                Console.WriteLine("\n");
                Console.WriteLine("Вектор нев'язок: ");
                Console.WriteLine(String.Join("\n", ls.UVector));
            }

            catch (Exception e)
            {
                Console.WriteLine(e.Message);
            }

        }
    }

    public class GaussSolutionNotFound : Exception
    {
        public GaussSolutionNotFound(string msg)
            : base("Рішення не може бути знайдено: \r\n" + msg)
        {
        }
    }

    public class LinearSystem
    {
        private double[,] initial_a_matrix;
        private double[,] a_matrix; // матриця A
        private double[] x_vector;  // вектор невідомих x
        private double[] initial_b_vector;
        private double[] b_vector;    // вектор b
        private double[] u_vector;    // вектор нев'язок U
        private double eps;           // порядок точності для дійсних чисел
        private int size;             // розмірність завдання


        public LinearSystem(double[,] a_matrix, double[] b_vector)
            : this(a_matrix, b_vector, 0.0001)
        {
        }
        public LinearSystem(double[,] a_matrix, double[] b_vector, double eps)
        {
            if (a_matrix == null || b_vector == null)
                throw new ArgumentNullException("Один з параметрів дорівнює null.");

            int b_length = b_vector.Length;
            int a_length = a_matrix.Length;
            if (a_length != b_length * b_length)
                throw new ArgumentException(@"Кількість рядків та стовпців в матриці А повинно співпадати з кількістю елементів у векторі B.");

            this.initial_a_matrix = a_matrix;  // запамятовуємо вихідну матрицю
            this.a_matrix = (double[,])a_matrix.Clone(); // з її копіює будемо проводити обчислення
            this.initial_b_vector = b_vector;  // запамятовуємо вихідний вектор
            this.b_vector = (double[])b_vector.Clone();  // з його копією будемо проводити обчислення
            this.x_vector = new double[b_length];
            this.u_vector = new double[b_length];
            this.size = b_length;
            this.eps = eps;

            GaussSolve();
        }

        public double[] XVector
        {
            get
            {
                return x_vector;
            }
        }

        public double[] UVector
        {
            get
            {
                return u_vector;
            }
        }

        // ініціалізація масива індексів стовпців
        private int[] InitIndex()
        {
            int[] index = new int[size];
            for (int i = 0; i < index.Length; ++i)
                index[i] = i;
            return index;
        }

        // пошук головного елементу в рядку матриці
        private double FindR(int row, int[] index)
        {
            int max_index = row;
            double max = a_matrix[row, index[max_index]];
            double max_abs = Math.Abs(max);
            //if(row < size - 1)
            for (int cur_index = row + 1; cur_index < size; ++cur_index)
            {
                double cur = a_matrix[row, index[cur_index]];
                double cur_abs = Math.Abs(cur);
                if (cur_abs > max_abs)
                {
                    max_index = cur_index;
                    max = cur;
                    max_abs = cur_abs;

                }
            }

            if (max_abs < eps)
            {
                if (Math.Abs(b_vector[row]) > eps)
                    throw new GaussSolutionNotFound("Система рівнянь несумісна.");
                else
                    throw new GaussSolutionNotFound("Система рівнянь маєм множину розвязків.");
            }

            // зміна місцями індексів стовпців
            int temp = index[row];
            index[row] = index[max_index];
            index[max_index] = temp;

            return max;
        }

        // знаходження розвязку СЛАР методом Гаусса 
        private void GaussSolve()
        {
            int[] index = InitIndex();
            PrintMatrix(initial_a_matrix);
            GaussForwardStroke(index);
            GaussBackwardStroke(index);
            GaussDiscrepancy();
        }

        // Вивід початкової матриці
        private void PrintMatrix(double[,] initial_a_matrix)
        {
            Console.ForegroundColor = ConsoleColor.Green;
            Console.WriteLine("Початкова матриця А та вектор В");
            Console.WriteLine();
            for (int i = 0; i < initial_a_matrix.GetLength(0); ++i)
            {
                for (int j = 0; j < initial_a_matrix.GetLength(1); ++j)
                {
                    Console.Write(initial_a_matrix[i, j] + "\t");
                }
                Console.Write(initial_b_vector[i] + "\t");
                Console.WriteLine(); // Перехід на нову строку.
            }
            Console.WriteLine();
        }
        // Прямий хід методу Гаусса
        private void GaussForwardStroke(int[] index)
        {
            // переміщення по кожному рядку з верху до низу
            for (int i = 0; i < size; ++i)
            {
                // 1) вибір головного елементу
                double r = FindR(i, index);

                // 2) перетворення поточного рядку матриці А
                for (int j = 0; j < size; ++j)
                    a_matrix[i, j] /= r;

                // 3) перетворення i-го елементу вектора b
                b_vector[i] /= r;

                // 4) Віднімання поточного рядка з усіх рядків, що під ним. 
                for (int k = i + 1; k < size; ++k)
                {
                    double p = a_matrix[k, index[i]];
                    for (int j = i; j < size; ++j)
                        a_matrix[k, index[j]] -= a_matrix[i, index[j]] * p;
                    b_vector[k] -= b_vector[i] * p;
                    a_matrix[k, index[i]] = 0.0;

                }
                Console.WriteLine("Максимальний елемент по модулю в рядку  " + (i + 1) + "  {0:f2}\t", r);
                Console.WriteLine();
                for (int a = 0; a < size; a++)
                {
                    for (int b = 0; b < size; b++)
                    {
                        Console.Write("{0:f2}\t", a_matrix[a, b]);
                    }
                    Console.Write("{0:f2}\t", b_vector[a]);
                    Console.WriteLine();
                }
                Console.WriteLine();
            }
        }

        // Зворотній хід метода Гаусса
        private void GaussBackwardStroke(int[] index)
        {
            // переміщення по кожному рядку знизу вверх
            for (int i = size - 1; i >= 0; --i)
            {
                // 1) задається початкове значення елементу x
                double x_i = b_vector[i];

                // 2) корегування даного значення
                for (int j = i + 1; j < size; ++j)
                    x_i -= x_vector[index[j]] * a_matrix[i, index[j]];
                x_vector[index[i]] = x_i;
            }
        }

        // Обчислення вектора невязок
        // U = b - x * A
        // x - розязок рівнняння, отриманого методом Гаусса
        private void GaussDiscrepancy()
        {
            for (int i = 0; i < size; ++i)
            {
                double actual_b_i = 0.0;   // результат перемноження i-строки 
                // вихідної матриці на вектор x
                for (int j = 0; j < size; ++j)
                    actual_b_i += initial_a_matrix[i, j] * x_vector[j];
                // i-й елемент вектора невязок
                u_vector[i] = initial_b_vector[i] - actual_b_i;
            }
        }

    }
}