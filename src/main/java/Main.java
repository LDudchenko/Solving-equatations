

public class Main {
    public static void main(String[] args) {
        Equation ob = new Equation();
        double apacheValue;
        ob.buildGraph1();
        ob.buildGraph2();

        System.out.println("Задание №1(y=x^2+sin(x)):");
        System.out.println("Первый корень:");

        ob.halvingMethod(-0.2, 0.3);
        ob.printAnswer("метода половинного деления");
        apacheValue = ob.bisectionMethod(-0.2, 0.3);
        printApacheValues(apacheValue);
        ob.compareValues(apacheValue);
        System.out.println();

        ob.chordMethod(-0.2, 0.3);
        ob.printAnswer("метода хорд");
        apacheValue = ob.secantMethod(-0.2, 0.3);
        printApacheValues(apacheValue);
        ob.compareValues(apacheValue);
        System.out.println();

        ob.tangentMethod(-0.2, 0.3);
        ob.printAnswer("метода касательных");
        apacheValue = ob.newtonRaphson(-0.2, 0.3);
        printApacheValues(apacheValue);
        ob.compareValues(apacheValue);
        System.out.println();

        System.out.println("Второй корень:");
        ob.halvingMethod(-1.25, -0.5);
        ob.printAnswer("метода половинного деления");
        apacheValue = ob.bisectionMethod(-1.25, -0.5);
        printApacheValues(apacheValue);
        ob.compareValues(apacheValue);
        System.out.println();

        ob.chordMethod(-1.25, -0.5);
        ob.printAnswer("метода хорд");
        apacheValue = ob.secantMethod(-1.25, -0.5);
        printApacheValues(apacheValue);
        ob.compareValues(apacheValue);
        System.out.println();

        ob.tangentMethod(-1.25, -0.5);
        ob.printAnswer("метода касательных");
        apacheValue = ob.newtonRaphson(-1.25, -0.5);
        printApacheValues(apacheValue);
        ob.compareValues(apacheValue);
        System.out.println();
        System.out.println();

        System.out.println("Задание №2(y=(x-40)-sin(x-40)):");
        ob.iterationMethod(38, 42.2);
        ob.printAnswer("метода итераций");
        ob.compareValues(39.72);
        System.out.println("\n");

        SystemOfEquatation ob1 = new SystemOfEquatation();
        double A[][] = new double[][]{{5, 6, 7, 0}, {10, 10, -1, 1}, {15, 4, -3, 2}};

        ob1.buildGraph3();

        System.out.println("Задание №3:");
        ob1.Gauss(A);
        Root[] rootsGauss = ob1.getGaussRoots();
        for(int i=0; i< rootsGauss.length;i++){
            rootsGauss[i].print();
        }
        System.out.println();
        double[] values = new double[]{(double) 53/405, (double)-1/27, (double)-5/81};
        ob1.compareGaussValues(values);
        System.out.println("\n");

        System.out.println("Задание №4:");
        double[] rootsZeydel = ob1.Zeydel(new double[]{15.8,0.2});
        ob1.print(rootsZeydel, "Метод Зейделя: ");
        ob1.compareValues(rootsZeydel, new double[]{16, 0});
        System.out.println();

        double[] rootsNewton = ob1.Newton(new double[]{15.8, 0.1});
        ob1.print(rootsNewton, "Метод Ньютона: ");
        ob1.compareValues(rootsNewton, new double[]{16, 0});
    }

    public static void printApacheValues(double value) {
        System.out.printf("Значение полученное с помощью математической библиотеки %.4f.\n", value);
    }


}


