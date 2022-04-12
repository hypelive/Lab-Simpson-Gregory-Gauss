#include <iostream>
#include <cmath>
#include <vector>

double f(double x)
{
    return cos(x + x * x);
}

double SimpsonSum(double h, int m, std::vector<double>& functionValues)
{
    double sum = 0.0;
    for (int i = 0; i < m; i++)
    {
        int arrayIndex = 2 * i;
        sum += h / 6 * (functionValues[arrayIndex] + 4 * functionValues[arrayIndex + 1] + functionValues[arrayIndex + 2]);
    }

    return sum;
}

double SimpsonRemain(double h)
{
    return 2977.0/2880.0 * h * h * h * h;
}

void SimpsonResults()
{
    double a = 2;
    double b = 3;
    double h = 0.1;

    double SimpsonApproximations[]{ 0, 0, 0 };
    double SimpsonApproximationRemains[]{ 0, 0, 0 };
    double SimpsonRungeRemains[]{ 0, 0, 0 };
    double SimpsonRungeApproximations[]{ 0, 0, 0 };

    for (int i = 0; i < 3; i++)
    {
        int m = (b - a) / h;
        std::vector<double> functionValues;
        functionValues.resize(2 * m + 1);
        for (int i = 0; i < m; i++)
        {
            int arrayIndex = 2 * i;
            functionValues[arrayIndex] = f(a + h * i);
            functionValues[arrayIndex + 1] = f(a + h * (i + 0.5));
        }

        functionValues[2 * m] = f(b);


        SimpsonApproximations[i] = SimpsonSum(h, m, functionValues);
        SimpsonApproximationRemains[i] = SimpsonRemain(h);

        h *= 0.5;
    }

    SimpsonRungeRemains[0] = NAN;
    SimpsonRungeRemains[1] = 1.0 / 15.0 * (SimpsonApproximations[1] - SimpsonApproximations[0]);
    SimpsonRungeRemains[2] = 1.0 / 15.0 * (SimpsonApproximations[2] - SimpsonApproximations[1]);

    SimpsonRungeApproximations[0] = NAN;
    SimpsonRungeApproximations[1] = SimpsonApproximations[1] + SimpsonRungeRemains[1];
    SimpsonRungeApproximations[2] = SimpsonApproximations[2] + SimpsonRungeRemains[2];

    std::cout << "Simpson h results:\nSum: " << SimpsonApproximations[0] << "\nRem: " << SimpsonApproximationRemains[0] << "\nRungeRem:" << SimpsonRungeRemains[0] << "\nRungeSum:" << SimpsonRungeApproximations[0] << std::endl << std::endl;
    std::cout << "Simpson h/2 results:\nSum: " << SimpsonApproximations[1] << "\nRem: " << SimpsonApproximationRemains[1] << "\nRungeRem:" << SimpsonRungeRemains[1] << "\nRungeSum:" << SimpsonRungeApproximations[1] << std::endl << std::endl;
    std::cout << "Simpson h/4 results:\nSum: " << SimpsonApproximations[2] << "\nRem: " << SimpsonApproximationRemains[2] << "\nRungeRem:" << SimpsonRungeRemains[2] << "\nRungeSum:" << SimpsonRungeApproximations[2] << std::endl << std::endl;
}

// Euler-Maclaurin: http://neo-chaos.narod.ru/useful/nummethod/kalitkin_04_085-125.pdf
double GregorySum(double h, int m, std::vector<double>& functionValues)
{
    double sum = 0.0;
    for (int i = 0; i < m; i++)
    {
        sum += h / 2 * (functionValues[i] + functionValues[i + 1]);
    }
    sum += h / 24 * (-3 * functionValues[0] + 4 * functionValues[1] - functionValues[2] - functionValues[m - 2] + 4 * functionValues[m - 1] - 3 * functionValues[m]);

    return sum;
}

double GregoryRemain(double h)
{
    return 682.7 / 72.0 * h * h * h * h;
}

void GregoryResults()
{
    double a = 2;
    double b = 3;
    double h = 0.1;

    double GregoryApproximations[]{ 0, 0, 0 };
    double GregoryApproximationRemains[]{ 0, 0, 0 };
    double GregoryRungeRemains[]{ 0, 0, 0 };
    double GregoryRungeApproximations[]{ 0, 0, 0 };

    for (int i = 0; i < 3; i++)
    {
        int m = (b - a) / h;
        std::vector<double> functionValues;
        functionValues.resize(m + 1);
        for (int i = 0; i <= m; i++)
        {
            functionValues[i] = f(a + h * i);
        }

        GregoryApproximations[i] = GregorySum(h, m, functionValues);
        GregoryApproximationRemains[i] = GregoryRemain(h);

        h *= 0.5;
    }

    GregoryRungeRemains[0] = NAN;
    GregoryRungeRemains[1] = 1.0 / 15.0 * (GregoryApproximations[1] - GregoryApproximations[0]);
    GregoryRungeRemains[2] = 1.0 / 15.0 * (GregoryApproximations[2] - GregoryApproximations[1]);

    GregoryRungeApproximations[0] = NAN;
    GregoryRungeApproximations[1] = GregoryApproximations[1] + GregoryRungeRemains[1];
    GregoryRungeApproximations[2] = GregoryApproximations[2] + GregoryRungeRemains[2];

    std::cout << "Gregory h results:\nSum: " << GregoryApproximations[0] << "\nRem: " << GregoryApproximationRemains[0] << "\nRungeRem:" << GregoryRungeRemains[0] << "\nRungeSum:" << GregoryRungeApproximations[0] << std::endl << std::endl;
    std::cout << "Gregory h/2 results:\nSum: " << GregoryApproximations[1] << "\nRem: " << GregoryApproximationRemains[1] << "\nRungeRem:" << GregoryRungeRemains[1] << "\nRungeSum:" << GregoryRungeApproximations[1] << std::endl << std::endl;
    std::cout << "Gregory h/4 results:\nSum: " << GregoryApproximations[2] << "\nRem: " << GregoryApproximationRemains[2] << "\nRungeRem:" << GregoryRungeRemains[2] << "\nRungeSum:" << GregoryRungeApproximations[2] << std::endl << std::endl;
}

void GaussResults()
{
    double x0 = (5.0 - sqrt(1.0 / 3.0)) / 2.0;
    double x1 = (5.0 + sqrt(1.0 / 3.0)) / 2.0;

    double A1 = (2.5 - x0) / (x1 - x0);
    double A0 = 1.0 - A1;

    bool AST_not_less_2 = abs(0.5 * (x0 * x0 + x1 * x1) - 19.0 / 3.0) < 1e-8;
    bool AST_not_less_3 = abs(0.5 * (x0 * x0 * x0 + x1 * x1 * x1) - (81.0 - 16.0) / 4.0) < 1e-8;

    double approx = A0 * f(x0) + A1 * f(x1);
    double remainIntergal = 0.2 * (243.0 - 32.0) - 0.5 * (x1 + x0)*(81 - 16) + 1.0 / 3.0 * (x1 * x1 + 4 * x0 * x1 + x0 * x0) * 19 - x0 * x1 * (x0 + x1) * 5 + x0 * x0 * x1 * x1;

    std::cout << "Gauss results:\nPoints: " << x0 << " " << x1 << "\nCoefs: " << A0 << " " << A1 << "\nAST 2 check:" << AST_not_less_2 << "\nAST 3 check:" << AST_not_less_3 << "\nApprox: " << approx << "\nIntegral for remain: " << remainIntergal << std::endl;
}

int main()
{
    SimpsonResults();
    GregoryResults();
    GaussResults();
}
