import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
import org.apache.commons.math3.analysis.differentiation.UnivariateDifferentiableFunction;
import org.apache.commons.math3.analysis.solvers.BisectionSolver;
import org.apache.commons.math3.analysis.solvers.NewtonRaphsonSolver;
import org.apache.commons.math3.analysis.solvers.SecantSolver;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.knowm.xchart.*;
import org.knowm.xchart.style.Styler;
import org.knowm.xchart.style.markers.SeriesMarkers;

public class Equation {
    double epsilon = 0.01;
    double answer;
    int g = 4;
    int k = 0;

    public void halvingMethod(double a, double b) {
        if (solveEquation(a) * solveEquation(b) >= 0) {
            System.out.println("На выбранном промежутке нет корней!");
            return;
        }
        double c = (a + b) / 2;
        boolean continuation = true;
        do {
            if ((b - a) <= epsilon) {
                answer = c;
                continuation = false;
            }
            if (solveEquation(c) == 0) {
                answer = c;
                continuation = false;
            } else if (solveEquation(a) * solveEquation(c) < 0) {
                b = c;
            } else if (solveEquation(b) * solveEquation(c) < 0) {
                a = c;
            }
            c = (a + b) / 2;
        } while (continuation);

    }

    public void chordMethod(double a, double b) {
        if (solveEquation(a) * solveEquation(b) >= 0) {
            System.out.println("На выбранном промежутке нет корней!");
            return;
        }
        double xt = a - (solveEquation(a) * (b - a)) / (solveEquation(b) - solveEquation(a));
        double xp;
        boolean continuation = true;
        do {
            if (solveEquation(xt) == 0) {
                answer = xt;
                return;
            } else if (solveEquation(a) * solveEquation(xt) < 0) {
                b = xt;
            } else if (solveEquation(b) * solveEquation(xt) < 0) {
                a = xt;
            }
            xp = xt;
            xt = a - (solveEquation(a) * (b - a)) / (solveEquation(b) - solveEquation(a));
            if (Math.abs(xt - xp) < epsilon) {
                answer = xt;
                continuation = false;
            }
        } while (continuation);
    }

    public void tangentMethod(double a, double b) {
        if (solveEquation(a) * solveEquation(b) >= 0) {
            System.out.println("На выбранном промежутке нет корней!");
            return;
        }
        boolean continuation = true;
        double xt = 0, xp;
        if (solveEquation(a) * findDerivatives(a, 2) > 0) {
            xt = a;
        } else if (solveEquation(b) * findDerivatives(b, 2) > 0) {
            xt = b;
        }
        do {
            xp = xt;
            xt = xp - solveEquation(xp) / findDerivatives(xp, 1);
            if (Math.abs(xt - xp) < epsilon) {
                answer = xt;
                continuation = false;
            }
        } while (continuation);
    }

    public void iterationMethod(double a, double b){
        double xp, xt=(a+b)/2;
        boolean continuation = true;
        do{
            xp= xt;
            xt = solveEquation2(xp);
            if(Math.abs(xt-xp)<epsilon){
                answer = xt;
                continuation = false;
            }
        }while(continuation);
    }


    public double solveEquation(double x) {
        return Math.pow(x - g * k, 2) + Math.sin(x - g * k);
    }

    public double solveEquation2(double x) {
        return 40 + Math.sin(x-40);
        //return  0.2*Math.sin(x-2)+2;
    }

    public DerivativeStructure solveEquation(DerivativeStructure x) {

        return (x.subtract(g * k)).pow(2).add((x.subtract(g * k)).sin());
    }

    public void printAnswer(String str) {
        System.out.printf("Корнем уравнения, найденного с помощью %s, является %.4f.\n", str, answer);
    }

    public double getAnswer() {
        return answer;
    }

    public double findDerivatives(double xRealValue, int order) {
        int params = 1;
        double res;
        DerivativeStructure x = new DerivativeStructure(params, order, 0, xRealValue);
        DerivativeStructure y = solveEquation(x);
        res = y.getPartialDerivative(order);
        return res;
    }

    public double bisectionMethod(double a, double b){
        BisectionSolver solver = new BisectionSolver();
        UnivariateFunction func = (double x) -> Math.pow(x,2)+Math.sin(x);
        return solver.solve(1000, func, a, b);
    }

    public double secantMethod(double a, double b){
        SecantSolver solver = new SecantSolver();
        UnivariateFunction func = (double x) -> Math.pow(x,2)+Math.sin(x);
        return solver.solve(1000, func, a, b);
    }

    public double newtonRaphson(double a, double b){
        NewtonRaphsonSolver solver = new NewtonRaphsonSolver();
        UnivariateDifferentiableFunction f = new UnivariateDifferentiableFunction() {

            public double value(double x) {
                return Math.pow(x,2)+Math.sin(x);
            }

            @Override
            public DerivativeStructure value(DerivativeStructure t) throws DimensionMismatchException {
                return (t.pow(2)).add(t.sin());
            }
        };

        return solver.solve(1000, f, a, b);
    }

    public void compareValues(double apacheValue){
        System.out.printf("Разница значений равна %.4f.\n",Math.abs(answer-apacheValue));
    }

    public void buildGraph1() {
        double[] xData = new double[40];
        double[] yData = new double[40];

        XYChart chart = new XYChartBuilder().width(800).height(600).title("Task1").build();

        chart.getStyler().setLegendPosition(Styler.LegendPosition.OutsideS);
        chart.getStyler().setLegendLayout(Styler.LegendLayout.Horizontal);
        chart.getStyler().setZoomEnabled(true);
        chart.getStyler().setCursorEnabled(true);

        double j = -2.5;
        for (int i = 0; i < 40; i++) {
            xData[i] = j;
            yData[i] = solveEquation(xData[i]);
            j += 0.125;
        }

        XYSeries series = chart.addSeries( "y(x)=x^2+sin(x)", xData, yData);
        series.setMarker(SeriesMarkers.NONE);

        new SwingWrapper(chart).displayChart();
    }

    public void buildGraph2() {
        double[] xData = new double[40];
        double[] yData = new double[40];

        XYChart chart = new XYChartBuilder().width(800).height(600).title("Task2").build();

        chart.getStyler().setLegendPosition(Styler.LegendPosition.OutsideS);
        chart.getStyler().setLegendLayout(Styler.LegendLayout.Horizontal);
        chart.getStyler().setZoomEnabled(true);
        chart.getStyler().setCursorEnabled(true);

        double j = 37;
        for (int i = 0; i < 40; i++) {
            xData[i] = j;
            yData[i] = (xData[i]-40)-Math.sin(xData[i]-40);
            j += 0.125;
        }

        XYSeries seriesHelpFunc = chart.addSeries("x=40+sin(x-40)",  xData, yData);
        seriesHelpFunc.setMarker(SeriesMarkers.NONE);


        new SwingWrapper(chart).displayChart();
    }
}



