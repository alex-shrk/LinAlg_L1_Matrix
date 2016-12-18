import Exceptions.MatrixSizeError;
import Matr.LUMatrix;
import Matr.Matrix;
import Matr.SparseMatrix;

import java.util.ArrayList;
import java.util.Random;
import java.util.Scanner;

public class Lab2 {
    public static double getRandomValue(){
        double rangeMin = -100d;
        double rangeMax = 100d;
        return rangeMin + (rangeMax - rangeMin) * Math.random();
    }


    public static void main(String[] args) {

        try {
            Random random = new Random();
            Scanner sc = new Scanner(System.in);
            System.out.println("Введите размерность матрицы: ");
            int size_matrix = sc.nextInt();

            SparseMatrix matrix = new SparseMatrix(size_matrix);

            for (int i = 0; i < size_matrix; i++) {
                double sum = 0.;
                double value = getRandomValue();
                sum += Math.abs(value);

                matrix.setElement(i, size_matrix - 1 - i, value);
                ArrayList<Integer> indexes = new ArrayList<>();
                int count = 1 + random.nextInt(10);

                while (indexes.size() < count) {
                    int ind = random.nextInt(size_matrix);
                    if (!indexes.contains(ind) && (ind != i)) {
                        indexes.add(ind);
                    }
                }

                for (int k = 0; k < count; k++) {
                    value = getRandomValue();
                    sum += Math.abs(value);
                    matrix.setElement(i, indexes.get(k), value);
                }
                matrix.setElement(i, i, (3229* sum + getRandomValue()));
            }

            matrix.outputMatrixPortraitToTxt("portrait");
            double condition = matrix.conditionMatrix();
            System.out.println("Condition=" + condition);
            int filling_start = matrix.fillingOfMatrix();
            System.out.println("Заполнение начальной матрицы = " + filling_start);

            Matrix x=new Matrix(size_matrix,1);
            for (int i=0;i<size_matrix;i++){
                x.setElement(i,0,(double)i+1.);
            }

            SparseMatrix b = matrix.multiply(x);

            System.out.println("Прямой метод (LU):");
            LUMatrix LU = new LUMatrix(matrix, true,filling_start);
            Matrix solveLU = LU.solveSystem(b);

            System.out.println("Вычисление относительной погрешности решения СЛАУ");
            double d1 = x.subtraction(solveLU).euclideanNorm();
            double d2 = x.euclideanNorm();
            double delta = d1 / d2;
            System.out.println("    delta=: " + delta);

            System.out.println("Введите точность для метода Зейделя");
            double precision_Seidel = sc.nextDouble();
            Matrix solveSeidel = matrix.solveSystemSeidelMethod(b, precision_Seidel);
            int filling_seidel = matrix.fillingOfMatrix();
            System.out.println("Заполнение матрицы для метода Зейделя= " + (filling_seidel-filling_start));

            System.out.println("Вычисление относительной погрешности решения СЛАУ");
            d1 = x.subtraction(solveSeidel).euclideanNorm();
            delta = d1 / d2;
            System.out.println("    delta=: " + delta);


        } catch (MatrixSizeError matrixSizeError) {
            matrixSizeError.printStackTrace();
        }
    }
}
