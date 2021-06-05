import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.knowm.xchart.SwingWrapper;
import org.knowm.xchart.XYChart;
import org.knowm.xchart.XYChartBuilder;
import org.knowm.xchart.XYSeries;
import org.knowm.xchart.style.Styler;
import org.knowm.xchart.style.markers.SeriesMarkers;

import java.util.Arrays;

public class SystemOfEquatation {
    double epsilon =0.01;
    Root[] Xarr;

    public void Gauss(double[][] A) {

        Xarr = new Root[A.length];
        for(int j=0; j< Xarr.length; j++){
            Xarr[j] = new Root(j+1,0);
        }

        for (int k = 0; k < A.length; k++) {
            double maxEl = Math.abs(A[k][k]);
            int iOfMaxEl = k, jOfMaxEl = k;
            for (int i = k; i < A.length; i++) {
                for (int j = k; j < A[0].length-1; j++) {
                    if (Math.abs(A[i][j]) > maxEl) {
                        maxEl = Math.abs(A[i][j]);
                        iOfMaxEl = i;
                        jOfMaxEl = j;
                    }
                }
            }

            if (iOfMaxEl != k) {
                double[] tempArr = A[k];
                A[k] = A[iOfMaxEl];
                A[iOfMaxEl] = tempArr;
            }
            if (jOfMaxEl != k) {
                for (int i = 0; i < A.length; i++) {
                    double temp = A[i][k];
                    A[i][k] = A[i][jOfMaxEl];
                    A[i][jOfMaxEl] = temp;
                }
                int temp = Xarr[jOfMaxEl].getNum();
                Xarr[jOfMaxEl].setNum(Xarr[k].getNum());
                Xarr[k].setNum(temp);
            }

            /*for (int i = 0; i < A.length; i++) {
                System.out.println(Arrays.toString(A[i]));
            }*/

            //Проверка есть ли решения у системы
            if(A[k][k]==0){
                for(int l=k; l<A.length; l++){
                    if(A[l][A[k].length-1]!=0){
                        System.out.println("Система не має розв'язків!");
                        return;
                    }
                }
                System.out.println("Система має безліч розв'язків!");
                return;
            }

            double mainEl = A[k][k];
            if(mainEl!=0) {
                for (int i = k; i < A[0].length; i++) {
                    A[k][i] /= mainEl;
                }
            }

            double Atemp[][] = A.clone();

            for (int i = k+1; i < A.length; i++) {
                double mult = Atemp[i][k];
                for (int j = k; j < A[0].length; j++) {
                    A[i][j] = Atemp[i][j]-Atemp[k][j]*mult;
                }
            }
        }


        //Обратный ход
        double sum;
        for(int k=0; k<A.length; k++){
            sum =0;
            for(int i=0; i<k; i++){
                sum += A[A.length-1-k][A[k].length-1-i-1]*Xarr[Xarr.length-i-1].getValue();
            }
            Xarr[Xarr.length-1-k].setValue((A[A.length-1-k][A[k].length-1]-sum)/A[A.length-1-k][A[k].length-1-(k+1)]);
        }
    }

    public double[] Zeydel(double[] x){
        double[] x_prev;
        double[] x_current=Arrays.copyOf(x,x.length);
        double max;
        do{
            x_prev = Arrays.copyOf(x_current, x_current.length);
            x_current[0]=funcForFirstEquatation(x_prev[0], x_prev[1]);
            x_current[1]=funcForSecondEquatation(x_current[0],x_prev[1]);
            max = Math.abs(x_current[0] - x_prev[0]);
            for (int i = 1; i < x.length; i++) {
                if (Math.abs(x_current[i] - x_prev[i]) > max) {
                    max = Math.abs(x_current[i] - x_prev[i]);
                }
            }
            if (max < epsilon) {
                break;
            }
        }while(true);
        return x_current;
    }

    public double funcForFirstEquatation(double x,double y){
        return 16-(Math.sin(x+y-16))/10;
    }

    public double funcForSecondEquatation(double x,double y){
        return Math.sin(x+y-16)/50;
    }

    public double firstEquatation(double[] x){
        return x[0]-16+(Math.sin(x[0]+x[1]-16))/10;
    }

    public double secondEquatation(double[] x){
        return x[1]-Math.sin(x[0]+x[1]-16)/50;
    }

    public double[] Newton(double[] x){
        double max;
        double[] X_curr = Arrays.copyOf(x,x.length);
        double[] X_prev;
        double[][] W = new double[x.length][x.length];

        do {
            X_prev = Arrays.copyOf(X_curr, X_curr.length);
            for (int i = 0; i < W.length; i++) {
                for (int j = 0; j < W.length; j++) {
                    W[i][j] = findPartialDerivatives(X_prev[0], X_prev[1], i, j);
                }
            }

            double[] F = new double[W.length];
            F[0] = firstEquatation(X_prev);
            F[1] = secondEquatation(X_prev);

            RealMatrix realMatrix = MatrixUtils.createRealMatrix(W);
            RealMatrix invertedMatrix = new LUDecomposition(realMatrix).getSolver().getInverse();
            X_curr = (MatrixUtils.createRealVector(X_prev)).subtract(invertedMatrix.operate((MatrixUtils.createRealVector(F)))).toArray();

            max = Math.abs(X_curr[0] - X_prev[0]);
            for (int i = 1; i < X_curr.length; i++) {
                if (Math.abs(X_curr[i] - X_prev[i]) > max) {
                    max = Math.abs(X_curr[i] - X_prev[i]);
                }
            }

            if (max <= epsilon) {
                break;
            }
        }while(true);
        return X_curr;

    }

    public DerivativeStructure firstEquation(DerivativeStructure x, DerivativeStructure y) {
        return  (x.subtract(16)).add(((x.add(y).subtract(16).sin()).divide(10)));


    }

    public DerivativeStructure secondEquation(DerivativeStructure x, DerivativeStructure y) {
        return y.subtract(((x.add(y).subtract(16)).sin()).divide(50));
    }

    public double findPartialDerivatives(double xRealValue, double yRealValue, int whichFunc, int whichArg) {
        int params = 2;
        int order = 2;
        DerivativeStructure x = new DerivativeStructure(params, order, 0, xRealValue);
        DerivativeStructure y = new DerivativeStructure(params, order, 1, yRealValue);
        DerivativeStructure g =(whichFunc == 0)?firstEquation(x,y):secondEquation(x,y);
        return (whichArg == 0)?g.getPartialDerivative(1, 0):g.getPartialDerivative(0, 1);
    }

    public Root[] getGaussRoots(){
        return Xarr;
    }

    public void compareGaussValues(double[] values){
        for(int i=0; i<values.length; i++){
            System.out.printf("Разница значений в %d корне равна %.4f.\n",Xarr[i].getNum(),Math.abs(Xarr[i].getValue()-values[Xarr[i].getNum()-1]));
        }
    }

    public void compareValues(double[] calculatedValues, double[] values){
        for(int i=0; i<values.length; i++){
            System.out.printf("Разница значений в %d корне равна %.4f.\n",i+1,Math.abs(calculatedValues[i]-values[i]));
        }
    }

    public void print(double[] roots, String str){
        System.out.print(str);
        for(int i=0; i<roots.length; i++){
            System.out.printf("x%d: %.4f; ", i+1,roots[i]);
        }
        System.out.println();
    }

    public void buildGraph3() {
        double[] xData = new double[40];
        double[] yData1 = new double[40];
        double[] yData2 = new double[40];

        XYChart chart = new XYChartBuilder().width(800).height(600).title("Task4").build();

        chart.getStyler().setLegendPosition(Styler.LegendPosition.OutsideS);
        chart.getStyler().setLegendLayout(Styler.LegendLayout.Horizontal);
        chart.getStyler().setZoomEnabled(true);
        chart.getStyler().setCursorEnabled(true);

        double j = 15.8;
        for (int i = 0; i < 40; i++) {
            xData[i] = j;
            yData1[i] = Math.asin(160-10*xData[i])-xData[i]+16;
            yData2[i] = Math.sin(xData[i]+yData1[i]-16)/50;
            j += 0.1;
        }

        XYSeries seriesHelpFunc = chart.addSeries("y=arcsin(160-10x)-x+16",  xData, yData1);
        seriesHelpFunc.setMarker(SeriesMarkers.NONE);

        XYSeries series = chart.addSeries("y=sin(x+arcsin(160-10x)-x+16-16)/50", xData, yData2);
        series.setMarker(SeriesMarkers.NONE);

        new SwingWrapper(chart).displayChart();
    }



}

class Root{
    private int num;
    private double value;

    public Root(int n, double v){
        num=n;
        value=v;
    }

    public Root(){
        num = 0;
    }

    public void setValue(double v){
        value =v;
    }

    public void setNum(int n){
        num=n;
    }

    public int getNum(){
        return num;
    }

    public double getValue(){
        return value;
    }

    public void print(){
        System.out.printf("x%d: %.4f; ",num,value);
    }
}
