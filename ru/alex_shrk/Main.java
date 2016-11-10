package ru.alex_shrk;

import Exceptions.MatrixSizeError;
import Matr.LUMatrix;
import Matr.Matrix;

import java.io.IOException;

public class Main {

    public static void main(String[] args) throws IOException, ClassNotFoundException, MatrixSizeError, CloneNotSupportedException {
        /*Scanner scanner=new Scanner(System.in);
        System.out.println("Input size of sq matrix");
        System.out.println("size:");
        int sizeSq=scanner.nextInt();

        Matrix matrixSq = new Matrix(sizeSq);
        matrixSq.setRange(0.0,10.0);
        matrixSq.generateMatrix();
        matrixSq.outputToTxt("matrSq");
        Matrix matrixSqTr=matrixSq.transpose();
        matrixSqTr.outputToTxt("matrSqTr");
        matrixSq.getVectorB().outputToTxt("matrSqVectB");
        List<Matrix> list=matrixSq.decompLU();
        list.get(0).outputToTxt("l_mSq");
        list.get(1).outputToTxt("u_mSq");
        System.out.println("Determinant"+matrixSq.determinant());*/



       /* System.out.println("Input size of matrix");
        System.out.println("m:");
        int m=scanner.nextInt();
        Matr vectorM=new Matr(m,1);
        vectorM.generateMatrix();
        vectorM.outputToTxt("vectorM.txt");
        Matr vectMTr=vectorM.transpose();
        vectMTr.outputToTxt("vectorMTr.txt");


        System.out.println("Input size of matrix");
        System.out.println("n:");
        int n_2=scanner.nextInt();
        Matr vectorN=new Matr(1,n_2);
        vectorN.generateMatrix();
        vectorN.outputToTxt("vectorN.txt");
        Matr vectNTr=vectorN.transpose();
        vectNTr.outputToTxt("vectorNTr.txt");*/

        /*System.out.println("Input size of matrix");
        System.out.println("size:");
        int size=scanner.nextInt();
        Matrix matr=new Matrix(size);
        matr.getSLAU().outputToTxt("slau");*/
        //matrix2.inputFromBin("matr1.bin");
        //matrix2.outputToTxt("matr2.txt");
        //Matr slau=matrix1.getSLAU();
        //slau.outputToTxt("slau.txt");

        /*
        Matrix slau=m33.getVectorB();
        slau.outputToTxt("slau33");
        List<Matrix> list=m33.decompLU();
        list.get(0).outputToTxt("l_m33");
        list.get(1).outputToTxt("u_m33");
        System.out.println("Determinant"+m33.determinant());

        LUMatrix luMatrix=new LUMatrix(m33);
        luMatrix.outputToTxt("luMatr");
        luMatrix.getA().outputToTxt("luA");
        luMatrix.getU().outputToTxt("luU");
        luMatrix.getL().outputToTxt("luL");
        luMatrix.inverse().outputToTxt("luInv");*/

       /* Matrix a_slau = new Matrix(new double[][]{{1., 0., 2.}, {-2., 3., 4.}, {1., 0., 8.}});
        List<Matrix> a_lu=a_slau.decompLU();
        a_lu.get(0).outputToTxt("slau_l");
        a_lu.get(1).outputToTxt("slau_u");
        List<Matrix> a_slau_solved=a_slau.solveSystem(a_lu);
        a_slau_solved.get(0).outputToTxt("slau_x");
        a_slau_solved.get(1).outputToTxt("slau_y");*/
       /*Matrix m1=new Matrix(new double[][]{{2., 0., 4., -1.}, {1., -1., 1.,0.}});
       Matrix m2=new Matrix(new double[][]{{2.}, {1.}, {0.}, {-2.}});
        Matrix.multiply(m1,m2).outputToTxt("res.txt");*/
        /*Matrix m33 = new Matrix(new double[][]{{2., 1., 1.}, {4., 3., 3.}, {8., 7., 4.}});
        LUMatrix lu=m33.decompLU();
        LUMatrix luLead=m33.decompLULead();
        lu.getLMatrix().outputToTxt("l_");
        lu.getUMatrix().outputToTxt("u_");
        luLead.getLMatrix().outputToTxt("l_lead_");
        luLead.getUMatrix().outputToTxt("u_lead_");

        m33.solveSystem(lu).outputToTxt("solve_lu");
        m33.solveSystem(luLead).outputToTxt("solve_lu_lead");*/
        Matrix m =new Matrix(4);
        m.setEps(16);
        m.generateMatrix();
        LUMatrix luMatrix=m.decompLU();

        luMatrix.outputToTxt("lu");
        System.out.println("Определитель матрицы:"+luMatrix.getDeterminant());
        m.solveSystem(luMatrix).outputToTxt("solve");
        System.out.println("Обусловленность матрицы: "+m.conditionMatrix(luMatrix));


        /*Matrix m = new Matrix(new double[][]{{2., 3., 3., 1.}, {3., 5., 3., 2.}, {5., 7., 6.,2.},{4.,4.,3.,1.}});

        Matrix b=new Matrix(new double[][]{{1.}, {2.}, {3.},{4.}});
        Matrix.solveSystem(m,0.0001).outputToTxt("solve");*/





    }
}
