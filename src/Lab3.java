import Matr.LUMatrix;
import Matr.Matrix;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.awt.*;
import java.util.Arrays;
import java.util.Scanner;

import static java.lang.Math.*;


public class Lab3 {
    private static double a = 0.;
    private static double b = 1.;
    private static int N = 5;
    private static double h = (b-a)/N;

    public static double K(double s,double t){
         return tanh(s+t);
    }
    public static double f(double t){
        return 3-t;
    }
    public static double Ai(int i) throws Exception {
        if (i==0 || i==N)
            return h/2.;
        else if (i>0 && i<N)
            return h;
        else
            throw new Exception(i+" не принадлежит границам [1,"+N+"]");
    }
    public static double alpha(int k){
        return pow(10,-k);
    }

    public static double calculateBeta(Matrix A,int variant, double alpha) throws Exception {
        //b*
        if (variant==1){
            double[][] matr = new double[A.getCountRows()][A.getCountColumns()];
            for (int i = 0; i < A.getCountRows(); i++) {
                for (int j = 0; j < A.getCountColumns(); j++) {
                    matr[i][j] = A.getElement(i, j);
                }
            }
            SingularValueDecomposition svd = new SingularValueDecomposition(new Array2DRowRealMatrix(matr));

            double[] svdValues = svd.getSingularValues();
            Arrays.sort(svdValues);
            double min = svdValues[0];
            double max = svdValues[svdValues.length-1];
            System.out.println("min sigma^2 = "+min);
            System.out.println("max sigma^2 = "+max);
            return sqrt(min/2. + alpha);
        }
        //b**
        else if (variant == 2)
            return sqrt(alpha);
        else
            throw new Exception("incorrect variant");
    }

    public static void main(String[] args) throws Exception {
        XYSeriesCollection xyDataSet = new XYSeriesCollection();
        int numberEquations = 0;
        Scanner sc=new Scanner(System.in);
        while (numberEquations==0 || numberEquations<(N+1)) {
            System.out.println("Введите число уравнений >="+(N+1));
            numberEquations=sc.nextInt();
        }



        Matrix x = new Matrix(numberEquations,1);
        for (int i=0;i<numberEquations;i++){
            x.setElement(i,0,a+i*h);
        }

        Matrix b = new Matrix(numberEquations,1);
        for (int i=0;i<numberEquations;i++){
            b.setElement(i,0,f(x.getElement(i)));
        }

        Matrix A = new Matrix(numberEquations,N+1);
        for (int i=0;i<numberEquations;i++){
            for (int j=0;j<N+1;j++) {
                A.setElement(i, j, Ai(j) * K ( x.getElement(i),x.getElement(j) ) );
            }
        }

        //System.out.println("Condition (A)=" +A.conditionMatrix());
        System.out.println();

        System.out.println("Метод нормальных уравнений\n");
        for (int k=0;k<3;k++){
            double alpha = alpha(k);
            Matrix A_regul = Matrix.regularization(A,alpha);
            Matrix L = Matrix.decompositionLLt(A_regul);
            Matrix phi_1 = Matrix.solveSystemLLt(L,A.transpose().multiply(b));

            XYSeries series = new XYSeries("alpha_"+k+" = " + alpha);
            for (int i = 0; i < N + 1; i++) {
                series.add(x.getElement(i), phi_1.getElement(i));
            }
            xyDataSet.addSeries(series);

            if (k==0)
                System.out.println("Condition (A_t A) = " + A.transpose().multiply(A).conditionMatrix()+"\n");
            System.out.println("alpha_" + k + " = " + alpha);

            System.out.println("Condition (A_t A + alpha"+"_"+k+" * E) = " + A_regul.conditionMatrix());

            //todo проверка cond(A)^2 = cond(A_t*A)

            System.out.println();
        }


        double alpha = alpha(0);
        System.out.println();
        System.out.println("Решение с помощью расширенной нормальной системы");
        System.out.println("alpha = " + alpha);
        System.out.println();

        for (int k = 1; k <= 2; k++) {
            double beta = calculateBeta(A, k, alpha);
            Matrix extendedA = Matrix.extensionMatrix(A, alpha, beta);
           // Matrix extendedB = new Matrix((N + 1) * 2, 1);
            Matrix extendedB = new Matrix(2*(A.getM()+A.getN() + 2), 1);
            extendedB.setValuesOfMatrix(0.);
            extendedB.insert(0, 0, b);
            Matrix phi_2 = new LUMatrix(extendedA).solveSystem(extendedB);//todo replace

            XYSeries series = new XYSeries("| alpha = " + alpha + ", beta_"+k+" = " + beta);
            for (int i = 0; i < N + 1; i++) {
                series.add(x.getElement(i), phi_2.getElement(N+1+i));
            }
            xyDataSet.addSeries(series);


            System.out.println("beta_" + (k+1) + " = " + beta);
            System.out.println("Condition of extended matrix A = " + extendedA.conditionMatrix());
            System.out.println();
        }

        JFreeChart chart1 = ChartFactory.createXYLineChart(
                "Графики приближенных значений phi, полученные с различными alpha, beta ", "X", "Phi",
                xyDataSet, PlotOrientation.VERTICAL, true, true, true
        );
        XYPlot plot=chart1.getXYPlot();
        plot.setBackgroundPaint(Color.DARK_GRAY);
        XYItemRenderer renderer = plot.getRenderer();
        renderer.setSeriesPaint(0, new Color(247, 255, 44));
        renderer.setSeriesPaint(1, new Color(255, 0, 255));
        renderer.setSeriesPaint(2, new Color(14, 147, 255));
        renderer.setSeriesPaint(3, new Color(14, 252, 255));
        renderer.setSeriesPaint(4, new Color(0, 255, 0));
        renderer.setStroke(new java.awt.BasicStroke(2));
        plot.setRenderer(renderer);

        JFrame frame1 = new JFrame("Graph");
        frame1.getContentPane().add(new ChartPanel(chart1));
        frame1.add(new ChartPanel(chart1));
        frame1.setSize(800, 600);
        frame1.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
        frame1.setVisible(true);

    }

}
