using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;

namespace MatrixMult
{
    class Matrix
    {
        int m; // количество строк
        public int Rows { get => m; }
        int n; // количество столбцов
        public int Columns {get => n;}
        double[][] data; // значения
        public double[][] AsDoubleArray { get => data; }
        public static Mutex mutex = new Mutex();
        public Matrix(int rows, int cols)
        {
            this.m = rows;
            this.n = cols;
            this.data = Create(m, n);
        }
        public Matrix(Matrix matrix)
        {
            this.m = matrix.Rows;
            this.n = matrix.Columns;
            this.data = Duplicate(matrix.AsDoubleArray);
        }
        public Matrix(double[][] matrix)
        {
            this.m = matrix.Length;
            this.n = matrix[0].Length;
            this.data = Duplicate(matrix);
        }
        static double[][] Create(int rows, int cols)
        {
            // Создаем матрицу, полностью инициализированную
            // значениями 0.0. Проверка входных параметров опущена.
            double[][] result = new double[rows][];
            for (int i = 0; i < rows; ++i)
                result[i] = new double[cols]; // автоинициализация в 0.0
            return result;
        }
        public void Randomize(double minVal, double maxVal, int seed)
        {
            // Возвращаем матрицу со значениями
            // в диапазоне от minVal до maxVal
            Random ran = new Random(seed);
            data = Create(Rows, Columns);
            for (int i = 0; i < Rows; ++i)
                for (int j = 0; j < Columns; ++j)
                    data[i][j] = (maxVal - minVal) * ran.NextDouble() + minVal;
        }
        public double Spur()
        {
            var size = Math.Min(Rows, Columns);
            double spur = 0;
            for (int i = 0; i < size; i++) spur += AsDoubleArray[i][i];
            return spur;
        }
        public bool AreEqual(Matrix matrix, double epsilon)
        {
            // True, если все значения в A равны
            // соответствующим значениям в B
            if (Columns != matrix.Columns && Rows != matrix.Rows)
                throw new Exception("Non-conformable matrices in MatrixAreEqual");
            for (int i = 0; i < Rows; ++i) // каждая строка A и B
                for (int j = 0; j < Columns; ++j) // каждый столбец A и B
                    if (Math.Abs(data[i][j] - matrix.AsDoubleArray[i][j]) > epsilon)
                        return false;
            return true;
        }
        public Matrix Product(Matrix matrix)
        {
            if (Columns != matrix.Rows)
                throw new Exception("Non-conformable matrices in MatrixProduct");
            double[][] result = Create(Rows, matrix.Columns);
            for (int i = 0; i < Rows; ++i) // каждая строка A
                for (int j = 0; j < matrix.Columns; ++j) // каждый столбец B
                    for (int k = 0; k < Columns; ++k)
                        result[i][j] += data[i][k] * matrix.AsDoubleArray[k][j];
            return new Matrix(result);
        }
        public Matrix ParallelProduct(Matrix matrix)
        {
            if (Columns != matrix.Rows)
                throw new Exception("Non-conformable matrices in MatrixProductParallel");
            double[][] result = Create(Rows, matrix.Columns);
            Parallel.For(0, Rows, i => {// каждая строка A
                for (int j = 0; j < matrix.Columns; ++j) // каждый столбец B
                    for (int k = 0; k < Columns; ++k)
                        result[i][j] += data[i][k] * matrix.AsDoubleArray[k][j];
            });
            return new Matrix(result);
        }
        private static void ProductPart(Matrix matrixA, Matrix matrixB, int start, int end, Matrix result)
        {
            for (int i = start; i < end; i++)
                for (int j = 0; j < matrixB.Columns; j++)
                    for (int k = 0; k < matrixA.Columns; k++)
                        result[i,j] += (matrixA[i,k] * matrixB[k,j]);
        }
        public Matrix ParallelProduct2(Matrix matrix, int threads)
        {
            if (Columns != matrix.Rows)
                throw new Exception("Non-conformable matrices in MatrixProductParallel2");
            Matrix result = new Matrix(Rows, matrix.Columns);
            List<Thread> ts = new List<Thread>();
            int start = 0;
            int len = Rows / threads + 1;
            for (int i = 0; i < threads; i++) 
            {
                var t = new Thread((object param) => {
                    var p = (object[])param;
                    ProductPart((Matrix)p[0], (Matrix)p[1], (int)p[2], (int)p[3], (Matrix)p[4]); 
                });
                ts.Add(t);
                t.Start(new object[5]{ this, matrix, start, (start + len > Rows) ? Rows : (start + len), result}) ;
                start+=len;
                //t.Join();
            }
            foreach (Thread t in ts) t.Join();
            return result;
        }
        public static Matrix operator *(Matrix a, Matrix b)
        {
            Matrix result = new Matrix(a.Rows, b.Columns);
            List<Thread> threads = new List<Thread>();
            for (int i = 0; i < a.Rows; i++)
                for (int j = 0; j < b.Columns; j++)
                {
                    int tempi = i;
                    int tempj = j;
                    Thread thread = new Thread(() => VectorMult(tempi, tempj, a, b, result));
                    thread.Start();
                    threads.Add(thread);
                }
            foreach (Thread t in threads)
                t.Join();
            return result;
        }
        public static void VectorMult(int tmpi, int tmpj, Matrix a, Matrix b, Matrix result)
        {
            mutex.WaitOne();
            int i = tmpi;
            int j = tmpj;
            double[] x = a.GetRow(i);
            double[] y = b.GetColumn(j);

            for (int k = 0; k < x.Length; k++)
                result[i, j] += x[k] * y[k];

            mutex.ReleaseMutex();
        }
        public double[] GetColumn(int i)
        {
            double[] res = new double[Rows];
            for (int j = 0; j < Rows; j++)
                res[j] = AsDoubleArray[j][i];
            return res;
        }
        public double[] GetRow(int i)
        {
            double[] res = new double[Columns];
            for (int j = 0; j < Columns; j++)
                res[j] = AsDoubleArray[i][j];
            return res;
        }
        public double this[int i, int j]
        {
            get { return AsDoubleArray[i][j]; }
            set { data[i][j] = value; }
        }
        static double[][] Duplicate(double[][] matrix)
        {
            // Предполагается, что матрица не нулевая
            double[][] result = Create(matrix.Length, matrix[0].Length);
            for (int i = 0; i < matrix.Length; ++i) // Копирование значений
                for (int j = 0; j < matrix[i].Length; ++j)
                    result[i][j] = matrix[i][j];
            return result;
        }
        public Matrix Decompose(out int[] perm, out int toggle)
        {
            // Разложение LUP Дулитла. Предполагается,
            // что матрица квадратная.
            double[][] result = Duplicate(AsDoubleArray);
            perm = new int[Rows];
            for (int i = 0; i < Rows; ++i) { perm[i] = i; }
            toggle = 1;
            for (int j = 0; j < Rows - 1; ++j) // каждый столбец
            {
                double colMax = Math.Abs(result[j][j]); // Наибольшее значение в столбце j
                int pRow = j;
                for (int i = j + 1; i < Rows; ++i)
                {
                    if (result[i][j] > colMax)
                    {
                        colMax = result[i][j];
                        pRow = i;
                    }
                }
                if (pRow != j) // перестановка строк
                {
                    double[] rowPtr = result[pRow];
                    result[pRow] = result[j];
                    result[j] = rowPtr;
                    int tmp = perm[pRow]; // Меняем информацию о перестановке
                    perm[pRow] = perm[j];
                    perm[j] = tmp;
                    toggle = -toggle; // переключатель перестановки строк
                }
                if (Math.Abs(result[j][j]) < 1.0E-20)
                    return null;
                for (int i = j + 1; i < Rows; ++i)
                {
                    result[i][j] /= result[j][j];
                    for (int k = j + 1; k < Rows; ++k)
                        result[i][k] -= result[i][j] * result[j][k];
                }
            } // основной цикл по столбцу j
            return new Matrix(result);
        }
        static double[] HelperSolve(double[][] luMatrix, double[] b)
        {
            // Решаем luMatrix * x = b
            int n = luMatrix.Length;
            double[] x = new double[n];
            b.CopyTo(x, 0);
            for (int i = 1; i < n; ++i)
            {
                double sum = x[i];
                for (int j = 0; j < i; ++j)
                    sum -= luMatrix[i][j] * x[j];
                x[i] = sum;
            }
            x[n - 1] /= luMatrix[n - 1][n - 1];
            for (int i = n - 2; i >= 0; --i)
            {
                double sum = x[i];
                for (int j = i + 1; j < n; ++j)
                    sum -= luMatrix[i][j] * x[j];
                x[i] = sum / luMatrix[i][i];
            }
            return x;
        }
        public Matrix Inverse()
        {
            double[][] result = Duplicate(AsDoubleArray);
            int[] perm;
            int toggle;
            double[][] lum = Decompose(out perm, out toggle).AsDoubleArray;
            if (lum == null)
                throw new Exception("Unable to compute inverse");
            double[] b = new double[Rows];
            for (int i = 0; i < Rows; ++i)
            {
                for (int j = 0; j < Rows; ++j)
                {
                    if (i == perm[j])
                        b[j] = 1.0;
                    else
                        b[j] = 0.0;
                }
                double[] x = HelperSolve(lum, b);
                for (int j = 0; j < Rows; ++j)
                    result[j][i] = x[j];
            }
            return new Matrix(result);
        }
        public double Determinant()
        {
            int[] perm;
            int toggle;
            Matrix lum = Decompose(out perm, out toggle);
            if (lum.AsDoubleArray == null)
                throw new Exception("Unable to compute MatrixDeterminant");
            double result = toggle;
            for (int i = 0; i < lum.Rows; ++i)
                result *= lum.AsDoubleArray[i][i];
            return result;
        }
        public string AsString()
        {
            string s = "";
            for (int i = 0; i < data.Length; ++i)
            {
                for (int j = 0; j < data[i].Length; ++j)
                    s += data[i][j].ToString("F3").PadLeft(8) + " ";
                s += Environment.NewLine;
            }
            return s;
        }
    }
}
