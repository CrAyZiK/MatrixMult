using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
namespace MatrixMult
{
    class Program
    {
        static void Main(string[] args)
        {
            int size = 400;
            Matrix matrixA = new Matrix(size, size);
            Matrix matrixB = new Matrix(size, size);
            matrixA.Randomize(-100, 100, (int)DateTime.UtcNow.Ticks);
            matrixB.Randomize(-100, 100, (int)DateTime.UtcNow.Ticks);
            var watch = System.Diagnostics.Stopwatch.StartNew();
            Matrix matrixC1 = matrixA.Product(matrixB);
            watch.Stop();
            Console.WriteLine($"Spur: {matrixC1.Spur()}");
            Console.WriteLine($"Time: {watch.ElapsedMilliseconds}ms");
            watch = System.Diagnostics.Stopwatch.StartNew();
            Matrix matrixC2 = matrixA.ParallelProduct2(matrixB, 12);
            watch.Stop();
            Console.WriteLine($"Spur: {matrixC2.Spur()}");
            Console.WriteLine($"Time: {watch.ElapsedMilliseconds}ms");
            double eps = 0.00001;
            Console.WriteLine($"MatrixC1 == MatrixC2 with eps ({eps}): {matrixC1.AreEqual(matrixC2, eps)}");
            //Console.WriteLine(matrixC1.AsString());
            //Console.WriteLine(matrixC2.AsString());
            Console.ReadKey();
        }
    }
}
