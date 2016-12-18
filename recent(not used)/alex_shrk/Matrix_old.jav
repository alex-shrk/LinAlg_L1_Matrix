package Matr;

import java.io.*;
import java.util.Random;

public class Matrix_old {
    private double[][] matrix;
    private int n;
    private String default_path="D://LinAlg//L1/matrix.txt";

    public Matrix_old() {
    }

    /** int n - size matrix */

    public Matrix_old(int n) {
        this.n=n;
        this.matrix = new double[n][n];
    }
    public void generateMatrix() {
        System.out.println("Generating matrix. Size="+n+"x"+n);
        Random random = new Random();
        for (int i = 0; i < n; i++) {
            for (int index = 0; index < n; index++) {
                matrix[i][index] = random.nextDouble();
            }
        }
        System.out.println("Complete");
    }
    public void inputFromFile(String name) throws IOException, ClassNotFoundException {
        FileInputStream fis=new FileInputStream(name);
        ObjectInputStream oin=new ObjectInputStream(fis);
        this.matrix= (double[][]) oin.readObject();
        this.n=matrix.length;

    }
    public void outputToBin(String name) throws IOException {
        FileOutputStream fos=new FileOutputStream(name);
        ObjectOutputStream oos= new ObjectOutputStream(fos);
        oos.writeObject(matrix);
        oos.flush();
        oos.close();
    }
    /*public void inputFromFile(String path) throws FileNotFoundException {
        System.out.println("Inputing from file...");
        File file = new File(path);
        if (!file.exists()){
            throw new FileNotFoundException(file.getName());
        }
        try{
            BufferedReader in = new BufferedReader(new FileReader(file.getAbsoluteFile()));
            String line;
            try{
                while ((line=in.readLine())!=null){
                    String[] value=line.split(" ");
                    for (int i = 0; i < n; i++) {
                        for (int index = 0; index < n; index++) {
                            matrix[i][index] = Double.parseDouble(value[i]);
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

    }
    */
    public  void outputToTxt(String path){
        System.out.println("Outputing to file...");
        if (path.isEmpty())
            path=this.default_path;
        File file = new File(path);
        try{
            if (!file.exists())
                file.createNewFile();

            PrintWriter out=new PrintWriter(file.getAbsoluteFile());
            try{
                for (int i = 0; i < n; i++) {
                    for (int index = 0; index < n; index++) {
                        out.print(matrix[i][index]);
                        out.print(' ');
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
