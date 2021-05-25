using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms.DataVisualization.Charting;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra.Double;

namespace InterpolationLab
{
    static class Interpolator
    {
        public const double a = 0;
        public const double b = 2;
        public const double drawA = a - 0;
        public const double drawB = b + 0;
        public const double delta = 0.01; //delta for drawing
        public const int N = 3; //nodes count
        private static double GetFunctionValue(double x)
        {
            return Math.Sin(x);
        }
        private static List<double> GetNodesCoords()
        {
            List<double> nodesCoords = new List<double>();
            nodesCoords.Add(a);
            for (int i = 1; i < N; i++)
            {
                nodesCoords.Add(a + (i * (b - a) / (N - 1)));
            }
            return nodesCoords;
        }
        private static List<double> GetXCoords()
        {
            List<double> xCoords = new List<double>();
            for (double i = drawA; i <= drawB; i += delta)
            {
                xCoords.Add(i);
            }
            return xCoords;
        }
        private static void GetGenFunctionXYValues(out List<double> xs, out List<double> ys)
        {
            xs = GetNodesCoords();
            ys = new List<double>();
            foreach (var x in xs)
            {
                ys.Add(GetFunctionValue(x));
            }
        }

        private static void GetFunctionXYValues(out List<double> xs, out List<double> ys, Polynomial f, double from, double to)
        {
            xs = new List<double>();
            ys = new List<double>();
            for (double i = from; i <= to; i += delta)
            {
                xs.Add(i);
                ys.Add(f.Evaluate(i));
            }
        }
        private static void DrawChart(string name, Chart chart, List<double> xs, List<double> ys)
        {
            chart.Series[name].Points.DataBindXY(xs, ys);
            chart.Invalidate();
        }
        public static void DrawGenFunction(Chart chart)
        {
            List<double> xs = new List<double>();
            List<double> ys = new List<double>();
            GetGenFunctionXYValues(out xs, out ys);
            DrawChart("Key", chart, xs, ys);
        }
        public static void DrawLagrangeFunction(Chart chart)
        {
            List<double> xs = new List<double>();
            List<double> ys = new List<double>();
            GetGenFunctionXYValues(out xs, out ys);
            GetFunctionXYValues(out xs, out ys, BuildLagrangePolynomial(xs, ys), drawA, drawB);
            DrawChart("Lagrange", chart, xs, ys);
        }
        public static void DrawNewtonFunction(Chart chart)
        {
            List<double> xs = new List<double>();
            List<double> ys = new List<double>();
            GetGenFunctionXYValues(out xs, out ys);
            GetFunctionXYValues(out xs, out ys, BuildNewtonPolynomial(xs, ys), drawA, drawB);
            DrawChart("Newton", chart, xs, ys);
        }
        public static void DrawStandartFunction(Chart chart)
        {
            List<double> xs = new List<double>();
            List<double> ys = new List<double>();
            GetGenFunctionXYValues(out xs, out ys);
            GetFunctionXYValues(out xs, out ys, BuildStandartPolynomial(xs, ys), drawA, drawB);
            DrawChart("Standart", chart, xs, ys);
        }
        public static void DrawSplineFunction(Chart chart)
        {
            List<double> xs = new List<double>();
            List<double> ys = new List<double>();
            GetGenFunctionXYValues(out xs, out ys);
            List<double> newX = null;
            List<double> newY = null;
            BuildSpline(xs, ys, out newX, out newY);
            DrawChart("Spline", chart, newX, newY);
        }
        private static Polynomial BuildLagrangePolynomial(List<double> xs, List<double> ys)
        {
            Polynomial l = new Polynomial(0, 0);
            for (int j = 0; j < xs.Count; j++)
            {
                Polynomial lb = BuildLagrangeBasicPolynomial(j, xs);
                l += ys[j] * lb;
            }
            return l;
        }
        private static Polynomial BuildLagrangeBasicPolynomial(int j, List<double> xs)
        {
            Polynomial lb = new Polynomial(1, 0);
            for (int m = 0; m < xs.Count; m++)
            {
                if (m != j)
                {
                    Polynomial numerator = new Polynomial(-xs[m], 1);
                    double denominator = xs[j] - xs[m];
                    lb *= numerator / denominator;
                }
            }
            return lb;
        }
        private static Polynomial BuildNewtonPolynomial(List<double> xs, List<double> ys)
        {
            Polynomial n = new Polynomial(0, 0);
            for (int j = 0; j < xs.Count; j++)
            {
                Polynomial nb = BuildNewtonBasicPolynomial(j, xs);
                double a = BuildNewtonBasicMultiplier(j, xs, ys);
                n += a * nb;
            }
            return n;
        }
        private static double BuildNewtonBasicMultiplier(int j, List<double> xs, List<double> ys)
        {
            double m = 0;
            if (j == 0)
            {
                m += ys[j];
            }
            else
            {
                for (int i = 0; i <= j; i++)
                {
                    double numerator = ys[i];
                    double denominator = 1;
                    for (int k = 0; k <= j; k++)
                    {
                        if (k != i)
                        {
                            denominator *= xs[i] - xs[k];
                        }
                    }
                    m += numerator / denominator;
                }
            }
            return m;
        }
        private static Polynomial BuildNewtonBasicPolynomial(int j, List<double> xs)
        {
            Polynomial nb = new Polynomial(1, 0);
            for (int i = 0; i <= j - 1; i++)
            {
                Polynomial numerator = new Polynomial(-xs[i], 1);
                nb *= numerator;
            }
            return nb;
        }
        private static Polynomial BuildStandartPolynomial(List<double> xs, List<double> ys)
        {
            Polynomial s = new Polynomial(SolveSystem(xs, ys));
            return s;
        }
        private static double[] SolveSystem(List<double> xs, List<double> ys)
        {
            int n = xs.Count;
            var m = Matrix.Build.Random(n, n);
            for (int i = 0; i < n; i++)
            {
                for (int j = n - 1; j >= 0; j--)
                {
                    m[i, j] = Math.Pow(xs[i], j);
                }
            }
            var y = Vector.Build.DenseOfEnumerable(ys);
            var x = m.Solve(y);

            double[] solution = new double[n];
            for (int i = n - 1; i >= 0; i--)
            {
                solution[i] = x[i];
            }
            return solution;
        }
        private static void BuildSpline(List<double> xs, List<double> ys, out List<double> resX, out List<double> resY)
        {
            int n = xs.Count;

            resX = new List<double>();
            resY = new List<double>();

            double[] h = new double[n - 1];

            for (int i = 1; i < h.Length; i++)
            {
                h[i] = xs[i] - xs[i - 1];
            }

            double[] a = new double[n - 1];

            for (int i = 0; i < a.Length; i++)
            {
                a[i] = ys[i];
            }

            double[] prC = SolveTridiagonal(h, a, n);
            double[] c = new double[n];
            c[0] = 0;
            for (int i = 1; i < prC.Length; i++)
            {
                c[i] = prC[i];
            }
            double[] d = new double[n];
            
            for(int i = 0; i < d.Length; i++)
            {
                d[i] = (c[i + 1] - c[i]) / (3 * h[i + 1]);
            }

            double[] b = new double[n - 1];

            for(int i = 0; i < b.Length; i++)
            {
                b[i] = (a[i + 1] - a[i]) / h[i + 1] + h[i + 1] * (2 * c[i + 1] + c[i]) / 3;
            }

            Polynomial[] polynomials = new Polynomial[n - 1];
            for(int i = 0; i < polynomials.Length; i++)
            {
                Polynomial xm = new Polynomial(-xs[i], 1);
                polynomials[i] = new Polynomial(0, 0);
                polynomials[i] += new Polynomial(a[i], 0);
                polynomials[i] += b[i] * xm;
                polynomials[i] += c[i] * xm * xm;
                polynomials[i] += d[i] * xm * xm * xm;

                for(double j = xs[i]; j < xs[i + 1]; j += delta)
                {
                    resX.Add(j);
                    resY.Add(polynomials[i].Evaluate(j));
                }
            }
        }

        private static double[] SolveTridiagonal(double[] h, double[] a, int n)
        {
            var m = Matrix.Build.Dense(n-2, n-2, 0.0);
            m[0, 0] = 2 * (h[0] + h[1]);
            m[0, 1] = h[1];
            m[n - 3, n - 3] = 2 * (h[n-3] + h[n-2]);
            m[n - 3, n - 4] = h[n - 2];
            for(int i = 1; i < m.ColumnCount-1; i++)
            {
                m[i, i] = 2 * (h[i] + h[i + 1]);
                m[i, i + 1] = h[i + 1];
                m[i, i - 1] = h[i];
            }
            var f = Vector.Build.Dense(n - 2);
            for (int i = 0; i < f.Count-1; i++)
            {
                f[i] = 3 * ((a[i + 2] - a[i + 1]) / h[i + 2] - (a[i + 1] - a[i]) / h[i + 1]);
            }
            var c = m.Solve(f).ToArray();
            return c;
        }

        public static void CubicInterpolation(List<double> sourceX, List<double> sourceY, List<double> newX, out List<double> newY)
        {
            int N = sourceX.Count;
            /*
             * Spline[i] = f[i] + b[i]*(x - x[i]) + c[i]*(x - x[i])^2 + d[i]*(x - x[i])^3
             * First: We prepare data for algorithm by calculate dx[i]. If dx[i] equal to zero then function return null.
             * Second: We need calculate coefficients b[i]. 
             * b[i] = 3 * ( (f[i] - f[i - 1])*dx[i]/dx[i - 1] + (f[i + 1] - f[i])*dx[i - 1]/dx[i] ),  i = 1, ... , N - 2
             * How calculate b[0] and b[N - 1] you can see below. And b can be find by means of tridiagonal matrix A[N, N].
             * 
             * A[N, N] - Tridiagonal Matrix:
             *      beta(0)     gama(0)     0            0           0   ...
             *      alfa(1)     beta(1)     gama(1)      0           0   ...
             *      0           alfa(2)     beta(2)     gama(2)      0
             *      ...
             * A*x=b
             * We calculate inverse of tridiagonal matrix by Gauss method and transforming equation A*x=b to the form I*x=b, where I - Identity matrix.
             * Fird: Now we can found coefficients c[i], d[i] where i = 0, ... , N - 2
             */

            int Nx = N - 1;
            double[] dx = new double[Nx];

            double[] b = new double[N];
            double[] alfa = new double[N];
            double[] beta = new double[N];
            double[] gama = new double[N];

            double[][] coefs = new double[4][];
            for (long i = 0; i < 4; i++)
                coefs[i] = new double[Nx];

            for (int i = 0; i + 1 <= Nx; i++)
            {
                dx[i] = sourceX[i + 1] - sourceX[i];
            }

            for (int i = 1; i + 1 <= Nx; i++)
            {
                b[i] = 3.0 * (dx[i] * ((sourceY[i] - sourceY[i - 1]) / dx[i - 1]) + dx[i - 1] * ((sourceY[i + 1] - sourceY[i]) / dx[i]));
            }

            b[0] = ((dx[0] + 2.0 * (sourceX[2] - sourceX[0])) * dx[1] * ((sourceY[1] - sourceY[0]) / dx[0]) +
                        Math.Pow(dx[0], 2.0) * ((sourceY[2] - sourceY[1]) / dx[1])) / (sourceX[2] - sourceX[0]);

            b[N - 1] = (Math.Pow(dx[Nx - 1], 2.0) * ((sourceY[N - 2] - sourceY[N - 3]) / dx[Nx - 2]) + (2.0 * (sourceX[N - 1] - sourceX[N - 3])
                + dx[Nx - 1]) * dx[Nx - 2] * ((sourceY[N - 1] - sourceY[N - 2]) / dx[Nx - 1])) / (sourceX[N - 1] - sourceX[N - 3]);

            beta[0] = dx[1];
            gama[0] = sourceX[2] - sourceX[0];
            beta[N - 1] = dx[Nx - 1];
            alfa[N - 1] = (sourceX[N - 1] - sourceX[N - 3]);
            for (long i = 1; i < N - 1; i++)
            {
                beta[i] = 2.0 * (dx[i] + dx[i - 1]);
                gama[i] = dx[i];
                alfa[i] = dx[i - 1];
            }
            double c = 0.0;
            for (long i = 0; i < N - 1; i++)
            {
                c = beta[i];
                b[i] /= c;
                beta[i] /= c;
                gama[i] /= c;

                c = alfa[i + 1];
                b[i + 1] -= c * b[i];
                alfa[i + 1] -= c * beta[i];
                beta[i + 1] -= c * gama[i];
            }

            b[N - 1] /= beta[N - 1];
            beta[N - 1] = 1.0;
            for (long i = N - 2; i >= 0; i--)
            {
                c = gama[i];
                b[i] -= c * b[i + 1];
                gama[i] -= c * beta[i];
            }

            for (int i = 0; i < Nx; i++)
            {
                double dzzdx = (sourceY[i + 1] - sourceY[i]) / Math.Pow(dx[i], 2.0) - b[i] / dx[i];
                double dzdxdx = b[i + 1] / dx[i] - (sourceY[i + 1] - sourceY[i]) / Math.Pow(dx[i], 2.0);
                coefs[0][i] = (dzdxdx - dzzdx) / dx[i];
                coefs[1][i] = (2.0 * dzzdx - dzdxdx);
                coefs[2][i] = b[i];
                coefs[3][i] = sourceY[i];
            }

            double[] Y = new double[newX.Count];
            int j = 0;
            for (int i = 0; i < N - 1; i++)
            {
                double h = 0.0;
                if (j >= newX.Count)
                    break;
                while (newX[j] < sourceX[i + 1])
                {
                    h = newX[j] - sourceX[i];
                    Y[j] = coefs[3][i] + h * (coefs[2][i] + h * (coefs[1][i] + h * coefs[0][i] / 3.0) / 2.0);
                    j++;
                    if (j >= newX.Count)
                        break;
                }
                if (j >= newX.Count)
                    break;
            }

            Y[Y.Length - 1] = sourceY[N - 1];
            newY = new List<double>();
            foreach(var e in Y)
            {
                newY.Add(e);
            }
        }
    }
}
