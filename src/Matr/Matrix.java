package Matr;

import Exceptions.MatrixSizeError;

import java.io.*;
import java.util.ArrayList;
import java.util.GregorianCalendar;
import java.util.List;
import java.util.Random;

public class Matrix {
    protected double[][] matrix;
    //m-count of rows
    //n-count of columns
    protected int m,n;
    protected double rangeMin=0.0;
    protected double rangeMax=10.0;
    protected double eps=0.0001;//todo сделать нормальное представление

    public Matrix() {
    }



    public int getM() {
        return m;
    }

    public int getN() {
        return n;
    }

    public void setRange(double rangeMin,double rangeMax){
        this.rangeMin=rangeMin;
        this.rangeMax=rangeMax;
    }
    public Matrix(int m, int n) {
        this.m=m;
        this.n=n;
        this.matrix = new double[m][n];
    }
    public Matrix(int sizeSquareMatr){
        this(sizeSquareMatr,sizeSquareMatr);
    }
    public Matrix(double[][] values){
        this(values.length,values[0].length);
        for (int i=0;i<m;i++){
            System.arraycopy(values[i],0,this.matrix[i],0,n);
        }
    }
    public void initializeValuesOfMatr(double val){
        for (int i=0;i<m;i++)
            for(int j=0;j<n;j++)
                matrix[i][j]=val;
    }

    public int getCountColumns() {
        return n;
    }

    public int getCountRows() {
        return m;
    }
    public double getElement(int m,int n) throws ArrayIndexOutOfBoundsException {
        try{
            return matrix[m][n];
        }
        catch (Exception e) {
            throw new ArrayIndexOutOfBoundsException("Element" + m + ' ' + n + "not found");
        }
    }
    public void setElement(int m,int n,double el) throws ArrayIndexOutOfBoundsException {
        try{
            matrix[m][n]=el;
        }
        catch (Exception e) {
            throw new ArrayIndexOutOfBoundsException("Element" + m + ' ' + n + "not found");
        }
    }
    public void setElement(int i, double el) throws ArrayIndexOutOfBoundsException {
        try{
            if (m==1)
                matrix[0][i]=el;
            else
            if (n==1)
                matrix[i][0]=el;
            else throw new MatrixSizeError("Невозможно присвоить элемент матрице. Только для векторов");
        }
        catch (Exception e) {
            throw new ArrayIndexOutOfBoundsException("Element" + m + ' ' + n + "not found");
        }
    }
    public Matrix multiply(Matrix matr) throws MatrixSizeError {
        return multiply(this,matr);
    }

    public static Matrix multiply(Matrix A, Matrix B) throws MatrixSizeError {
        int m_a=A.getCountRows();
        int n_a=A.getCountColumns();

        int m_b=B.getCountRows();
        int n_b=B.getCountColumns();


        if (n_a!=m_b)
            throw new MatrixSizeError("Произведение матриц невозможно: число столбцов 1-й матрицы!=числу строк 2-й матрицы");

        Matrix c=new Matrix(m_a,n_b);
        for(int i=0;i<c.m;i++){
            for(int j=0;j<c.n;j++){
                c.matrix[i][j]=0.0;
                for (int k=0;k<n_a;k++)
                    c.matrix[i][j]+=A.matrix[i][k]*B.matrix[k][j];
            }
        }
        return c;
    }
    public static Matrix addition(Matrix A, Matrix B) throws MatrixSizeError {
        int n_a=A.getCountColumns();
        int n_b=B.getCountColumns();

        int m_a=A.getCountRows();
        int m_b=B.getCountRows();

        if (m_a!=m_b && n_a!=n_b)
            throw new MatrixSizeError("");
        double el_a,el_b;

        Matrix C=new Matrix(m_a,n_a);


        for(int i=0;i<m_a;i++){
            for(int j=0;j<n_a;j++){
                el_a=A.getElement(i,j);
                el_b=B.getElement(i,j);
                C.setElement(i,j,el_a+el_b);
            }
        }
        return C;
    }

    public double getRandomDoubleValue(){
        Random random = new Random();
        return rangeMin+(rangeMax-rangeMin)*random.nextDouble();
    }
    public boolean isCorrectRange(double rangeMin,double rangeMax){
        return !Double.valueOf(rangeMax - rangeMin).isInfinite();
    }
    public boolean isCorrectRange(){
        return isCorrectRange(rangeMin,rangeMax);
    }



    public void generateMatrix() {
        System.out.println("Generating matrix. Size="+m+"x"+n);
        if (isCorrectRange()) {

            for (int i = 0; i < m; i++) {
                for (int j = 0; j < n; j++) {
                    matrix[i][j] = getRandomDoubleValue();
                }
            }
            System.out.println("Complete");
        }
        else throw new ArithmeticException("rangeMax-rangeMin=infinite");
    }

    public Matrix generateVectorX(){
        return generateVectorX(this.m);
    }
    public Matrix generateVectorX(int m){
        int n=1;
        Matrix x=new Matrix(m,n);
            for (int i=0;i<m;i++){
                x.setElement(i,i+1.);
            }

        return x;
    }

    public static Matrix transpose(Matrix matr){
        int m=matr.getCountRows();
        int n=matr.getCountColumns();
        Matrix transposedMatrix=new Matrix(matr.getCountColumns(),matr.getCountRows());
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                transposedMatrix.matrix[j][i] = matr.matrix[i][j];
            }
        }
        return transposedMatrix;

    }
    public Matrix transpose(){
        return transpose(this);
    }
    public void inputFromBin(String name) throws IOException, ClassNotFoundException {
        FileInputStream fis=new FileInputStream(name);
        ObjectInputStream oin=new ObjectInputStream(fis);
        this.matrix= (double[][]) oin.readObject();
        this.n=matrix.length;
        this.m=matrix[0].length;

    }
    public void outputToBin(String name) throws IOException {
        FileOutputStream fos=new FileOutputStream(name);
        ObjectOutputStream oos= new ObjectOutputStream(fos);
        oos.writeObject(matrix);
        oos.flush();
        oos.close();
    }

    /*public Matrix inverse() throws MatrixSizeError, MatrixNotExist {
        double det=determinant();
        if (det==0)
            throw new MatrixNotExist("Обратная матрица не существует");
        Matrix invMatr=new Matrix(m,n);
        for (int i=m-1;i>=0;i--){

        }


    }*/

    /** метод обрежет матрицу до размеров int order*/
    public static Matrix cropMatrix(int order,Matrix A){
        Matrix B=new Matrix(order);
        for (int i=0;i<order;i++){
            System.arraycopy(A.matrix[i],0,B.matrix[i],0,order);
        }
        return B;
    }

    public static double getMajorMinor(int order,Matrix A) throws Exception {
        //order-порядок главного минора матрицы А
        if (order>A.matrix.length)
            throw new Exception("порядок минора больше порядка матрицы");

        return determinant(cropMatrix(order, A));
    }

    public static boolean checkingPositivityOfMajorMinors(Matrix A) throws Exception {
        for (int i=1;i<A.matrix.length;i++) {
            if (getMajorMinor(i, A) <= 0) {
                return false;
            }
        }
        return true;

    }

    public static double determinant(Matrix A) throws MatrixSizeError {
        if (A.getM()!=A.getN())
            throw new MatrixSizeError("Матрица не квадратная");
        double det=1.0;
        List<Matrix> lu=A.decompLU();
        for (int i=0;i<lu.get(1).getCountColumns();i++)
            det*=lu.get(1).matrix[i][i];
        return det;
    }
    public static double determinant(List<Matrix> lu) throws MatrixSizeError {
        double det=1.0;
        for (int i=0;i<lu.get(1).getCountColumns();i++)
            det*=lu.get(1).matrix[i][i];
        return det;
    }


    public double determinant() throws MatrixSizeError {
        return determinant(this);
    }



    //Сферическая векторная норма

    /*public double norma2() {
        Matrix AA;
        if (m <= n) {
            AA = this.multiply(transpose());
        } else {
            AA = transpose().multiply(this);
        }

        Double max = 0.0;
        for (int k = 0; k < AA.m; k++) {
            max = Math.max(max, AA.data[k][k]);
        }
        return Math.sqrt(max);
    }*/

    public Matrix getVectorB() throws MatrixSizeError {
        return multiply(generateVectorX());
    }

    public List<Matrix> decompLU() throws MatrixSizeError {
        if (m!=n){
            throw new MatrixSizeError("Матрица не квадратная");
        }
        Matrix l=new Matrix(this.m,this.n);
        Matrix u=new Matrix(this.m,this.n);

        l.initializeValuesOfMatr(0.0);
        u.initializeValuesOfMatr(0.0);

        for (int i=0;i<m;i++){
            for (int j=i;j<n;j++){
                double val_u=matrix[i][j];
                for (int k=0;k<i;k++) {
                    val_u -= l.matrix[i][k] * u.matrix[k][j];
                }
                u.matrix[i][j]=val_u;

                double val_l=matrix[j][i];
                for (int k=0;k<i;k++){
                    val_l -= l.matrix[j][k] * u.matrix[k][i];
                }
                if (Math.abs(u.matrix[i][i])>=eps)
                    val_l/=u.matrix[i][i];
                else{
                    System.out.println("Элемент u"+i+" "+i+"меньше eps");
                    break;
                }


                l.matrix[j][i]=val_l;

            }
        }
        ArrayList<Matrix> lu=new ArrayList<>();
        lu.add(l);
        lu.add(u);
        return lu;

    }

    public List<Matrix> solveSystem(List<Matrix> lu){
        Matrix y=new Matrix(this.m,1);
        y.initializeValuesOfMatr(0.0);
        Matrix x=new Matrix(this.m,1);
        x.initializeValuesOfMatr(0.0);
        Matrix b = new Matrix(new double[][]{{6.,0.,0.}, {18.,0.,0.}, {24.,0.,0.}});
        //Matrix b=generateVectorX();

        int y_n=y.getN();
        int b_n= b.getN();
        double u=0.0;

        for (int i=0;i<=y_n;i++){
            for (int k=0;k<i;k++){
                y.matrix[i][0]-=lu.get(0).matrix[i][k]*y.matrix[k][0];
            }
            y.matrix[i][0]+=b.matrix[i][0];

            for(int k=i;k<y_n;k++){
                x.matrix[y_n-i][0]-=lu.get(1).matrix[y_n-i][k]*x.matrix[k][0];
            }
            u=lu.get(1).matrix[i][i];
            if (Math.abs(u)>=eps)
                x.matrix[y_n-i][0]+=y.matrix[i][0]/u;
            else{
                System.out.println("Элемент u"+i+" "+i+"меньше eps");
                break;
            }

        }
        List<Matrix> result=new ArrayList<>();
        result.add(x);
        result.add(y);
        return result;


    }

    /*public Matrix inverse() {
        double[][] result = new double[m][m];

        for (int p = m-1; p >= 0; p--) {
            Double val = 1.0;
            for (int k = p+1; k < m; k++) {
                val -= getElemU(p, k) * result[k][p];
            }
            val /= getElemU(p, p);
            result[p][p] = val;

            for (int i = p-1; i >= 0; i--) {
                val = 0.0;
                for (int k = i+1; k < m; k++) {
                    val -= getElemU(i, k) * result[k][p];
                }
                val /= getElemU(i, i);
                result[i][p] = val;
            }

            for (int j = p-1; j >= 0; j--) {
                val = 0.0;
                for (int k = j+1; k < m; k++) {
                    val -= result[p][k] * getElemL(k, j);
                }
                result[p][j] = val;
            }
        }

        return new Matrix(result);
    }*/





    /*public void inputFromFile(String path) throws FileNotFoundException {
            System.out.println("Inputing from file...");
            File file = new File(path);
            if (!file.exists()){
                throw new FileNotFoundException(file.getName());
            }
            try{
               // BufferedReader in = new BufferedReader(new FileReader(file.getAbsoluteFile()));
                Scanner sc=new Scanner(file);
                sc.useDelimiter(" ");
                try{
                    while (sc.hasNextDouble()){
                        String[] value=line.split(" ");
                        for (int i = 0; i < n; i++) {
                            for (int j = 0; j < n; j++) {
                                matrix[i][j] = Double.parseDouble(value[i]);
                            }
                        }
                    }
                } catch (IOException e) {
                    e.printStackTrace();
                }
                finally {
                    in.close();
                }
            }
            catch (IOException e){
                throw new RuntimeException(e);
            }
            System.out.println("Complete");

        }*/

    public  void outputToTxt(String path){
        System.out.println("Outputing to file...");
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
                        out.print(matrix[i][j]);
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
        System.out.println("Complete");
    }
}
