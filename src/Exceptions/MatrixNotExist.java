package Exceptions;

public class MatrixNotExist extends Exception {
    public MatrixNotExist() {
        super("Несовместимые размеры матриц");
    }
    public MatrixNotExist(String str){
        super(str);
    }
}
