package Matr;


import Exceptions.MatrixSizeError;

import java.io.*;
import java.util.GregorianCalendar;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.ListIterator;

public class SparseMatrix implements IMatrix {
    private int m;
    private int n;
    private LinkedList<Element>[] matrix;

    class Element{
        int index;
        double value;

        public Element(int index, double value) {
            this.index = index;
            this.value = value;
        }

        public int getIndex() {
            return index;
        }

        public double getValue() {
            return value;
        }
    }

    /** return vector (1,2,...m)*/
    public static SparseMatrix generateVector(int m){
        int n=1;
        SparseMatrix x=new SparseMatrix(m,n);
        for (int i=0;i<m;i++){
            x.setElement(i,0,i+1.);
        }

        return x;
    }


    public SparseMatrix(int m, int n) {
        this.m=m;
        this.n=n;
        matrix = new LinkedList[m];
        for (int i=0;i<m;i++){
            matrix[i]= new LinkedList<>();
        }

    }

    public SparseMatrix(int sizeSquareMatr) {
        this(sizeSquareMatr,sizeSquareMatr);
    }

    public SparseMatrix(double[][] values) {
        this(values.length,values[0].length);
        for (int i=0;i<m;i++){
            for(int j=0;j<n;j++)
                setElement(i,j,values[i][j]);
        }
    }

    @Override
    public int getM() {
        return m;
    }

    @Override
    public int getN() {
        return n;
    }

    

    @Override
    public int getCountColumns() {
        return n;
    }

    @Override
    public int getCountRows() {
        return m;
    }

    //return count of not null elements matrix
    public int fillingOfMatrix(){
        int count=0;
        for (int i=0;i<m;i++)
            count+=matrix[i].size();
        return count;
    }

    @Override
    public double getElement(int row, int column) {
        if (row<0 || row>=m)
            throw new IndexOutOfBoundsException("Incorrect bounds" + m + ' ' + "not in ["+0+","+(m-1)+"]");
        if (column<0 || column>=n)
            throw new IndexOutOfBoundsException("Incorrect bounds" + n + ' ' + "not in ["+0+","+(n-1)+"]");
        for (Element elem:matrix[row]){
            if (elem.index <column)
                continue;

            else if (elem.index >column) {
                break;
            }
            else if (elem.index ==column){
                return elem.value;
            }
        }
        return 0.d;
    }


    @Override
    public void setElement(int m, int n, double el) throws ArrayIndexOutOfBoundsException {
        if (m<0 || m>=this.m)
            throw new IndexOutOfBoundsException("Incorrect bounds" + m + ' ' + "not in ["+0+","+this.m+"]");
        if (n<0 || n>=this.n)
            throw new IndexOutOfBoundsException("Incorrect bounds" + n + ' ' + "not in ["+0+","+this.n+"]");

        ListIterator<Element> listIterator = matrix[m].listIterator();

        if (!matrix[m].isEmpty()){
            while(listIterator.hasNext()){
                Element elem=listIterator.next();
                if (elem.index <n)
                    continue;
                else if (elem.index >n) {
                    listIterator.previous();
                    break;
                }
                else if (elem.index ==n){
                    listIterator.remove();
                    break;
                }

            }
        }
        if (el!=0.d){
            Element valueElem = new Element(n,el);
            listIterator.add(valueElem);
        }

    }
    public SparseMatrix multiply(IMatrix B) throws MatrixSizeError {

        int m_b=B.getM();
        int n_b=B.getN();

        if (n!=m_b)
            throw new MatrixSizeError("Невозможно выполнить операцию умножения матриц: число столбцов 1-й матрицы!=числу строк 2-й матрицы");

        SparseMatrix c=new SparseMatrix(m,n_b);
        double value;
        for(int i=0;i<c.m;i++){
            for(int j=0;j<c.n;j++){
                value=0.d;
                for (Element elem:matrix[i])
                    value += elem.value * B.getElement(elem.index,j);
                if (value!=0.0d)
                    c.setElement(i,j,value);

            }
        }
        return c;
    }


    public double euclideanNorm() {
        double norm=0.0;
        for (int i=0;i<m;i++){
            for (Element elem:matrix[i]){
                norm+=Math.pow(Math.abs(elem.value), 2);
            }
        }
        return Math.sqrt(norm);
    }


    public double conditionMatrix() throws MatrixSizeError {
        LUMatrix luMatrix = new LUMatrix(this);
        double euclideanNormA=euclideanNorm();
        double euclideanNormInvA=luMatrix.getInverseMatrix().euclideanNorm();
        return euclideanNormA*euclideanNormInvA;
    }


    public Matrix solveSystemSeidelMethod(IMatrix b, double eps) throws MatrixSizeError {
        SparseMatrix LD = new SparseMatrix(m,n);//lower+diagonal matrix
        SparseMatrix R = new SparseMatrix(m,n);//upper matrix

        for (int i=0;i<m;i++) {
            for (Element el : matrix[i]) {
                if (i >= el.index) {
                    LD.setElement(i, el.index, el.value);
                } else {
                    R.setElement(i, el.index, el.value);
                }
            }
        }

        Matrix inverseLD = new LUMatrix(LD).getInverseMatrix();
        Matrix B = inverseLD.multiply(R).multiply(-1.);
        Matrix c = inverseLD.multiply(b);

        double q = B.euclideanNorm();


        Matrix x_k;
        Matrix x_k_1 = c;

        int count_of_iterations = 0;

        do{
            x_k = new Matrix(x_k_1);
            for (int i=0;i<m;i++){
                double value = b.getElement(i,0);
                for (int j=0;j<n;j++){
                    if (j!=i)
                        value-=getElement(i,j) * x_k_1.getElement(j);
                }
                double A_i_i = getElement(i,i);
                assert A_i_i!=0.d;

                x_k_1.setElement(i,value/A_i_i);
            }
            count_of_iterations++;
        }
        while (x_k_1.subtraction(x_k).euclideanNorm()> ((1-q)/q) * eps);

        System.out.println("Число итераций:"+count_of_iterations);
        return x_k_1;
    }



    public void outputMatrixPortraitToTxt(String path){
        int m=getM();
        int n=getN();

        StringBuilder portrait = new StringBuilder(m*(n+5));

        for (int i = 0; i < m; i++) {
            portrait.append("|||");
            LinkedList<Element> line = matrix[i];
            int k = 0;
            Iterator<Element> iterator = line.iterator();
            Element element = null;
            boolean skip = false;
            for (int j = 0; j < n; j++) {
                if (!skip && iterator.hasNext()) {
                    element = iterator.next();
                    skip = true;
                }
                if (element != null && element.index == j) {
                    portrait.append('+');
                    skip = false;
                } else {
                    portrait.append(' ');
                }
            }
            portrait.append("|||").append('\n');
        }

        try (FileOutputStream fileOutputStream = new FileOutputStream(path + ".txt");
             PrintStream file = new PrintStream(fileOutputStream)) {
            file.print(portrait);
            file.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }



    @Override
    public void outputToTxt(String path) {
        System.out.println("Writing to file...");
        String default_path="matrix"+ new GregorianCalendar().getTime().toString();
        if (path.isEmpty())
            path=default_path;
        File file = new File(path+".txt");
        try{
            if (!file.exists())
                file.createNewFile();

            PrintWriter out=new PrintWriter(file.getAbsoluteFile());
            try{

                for (int i = 0; i < m; i++) {
                    for (int j = 0; j < n; j++) {
                        out.print(getElement(i,j));
                        out.print('\t');
                    }
                    out.print('\n');
                }
            }
            finally {
                out.close();
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println("Complete\n");
    }

}
