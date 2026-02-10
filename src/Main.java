import java.util.*;

public class Main {
    public static void main(String[] args) {
        double[][] m = {{1,-3,-2,-4,0,0,0,0},{0,1,1,2,1,0,0,4},{0,2,0,3,0,1,0,5},{0,2,1,3,0,0,1,7}};
        double[] vals = {0,0,0,4,5,7};
        Simplex simplex = new Simplex(m, vals);
        double res = simplex.solve();
        System.out.println(res);
    }
}

class Simplex{
    Matrix matrix;
    double[] vals;
    Set<Integer> basics;
    Set<Integer> nonbasics;
    public Simplex(double[][] m, double[] vals){
        basics = new HashSet<>();
        nonbasics = new HashSet<>();
        matrix = new Matrix(m);
        this.vals = vals.clone();
        for(int i=0; i<vals.length; i++){
            if(vals[i]==0){
                nonbasics.add(i);
            }
            else{
                basics.add(i);
            }
        }
    }
    public double solve(){
        int pivot = choosePivot();
        while(pivot != -1){
            int kickRow = findKickRow(pivot);
            if(kickRow == -1){
                return Double.MIN_VALUE;
            }
            int[][] pivots = {{kickRow,pivot}};
            matrix.pivotify(pivots);
            pivot = choosePivot();
        }
        System.out.println(matrix);
        return matrix.vals[0][matrix.cols-1];
    }
    public int choosePivot(){
        for(int i=1; i<matrix.cols-1; i++){
            if(matrix.vals[0][i]<0){
                return i;
            }
        }
        return -1;
    }
    public int findKickRow(int pivot){
        double min = Double.MAX_VALUE;
        int kickRow = -1;
        for(int row=1; row<matrix.rows; row++){
            double ratio = matrix.vals[row][matrix.cols-1]/matrix.vals[row][pivot];
            if(ratio < 0) continue;
            if(ratio < min){
                min = ratio;
                kickRow = row;
            }
        }
        return kickRow;
    }
}

class Matrix{
    double[][] vals;
    int cols;
    int rows;
    public Matrix(double[][] v){
        vals = v.clone();
        cols = v[0].length;
        rows = v.length;
    }
    public Matrix(int a, int b){
        cols = a;
        rows = b;
        vals = new double[rows][cols];
    }
    public void pivotify(int[][] pivots){
        for(int[] pivot : pivots){
            clear(pivot);
            normalize(pivot);
        }
    }
    public void clear(int[] pivot){
        for(int row=0; row<rows; row++){
            if(pivot[0] == row) continue;
            double scalar = -vals[row][pivot[1]]/vals[pivot[0]][pivot[1]];
            rowAdd(scalar, pivot[0], row);
        }
    }
    public void normalize(int[] pivot){
        rowScale(1/vals[pivot[0]][pivot[1]], pivot[0]);
    }
    public void rowAdd(double scalar, int src, int dest){
        for(int col=0; col<cols; col++){
            vals[dest][col] += scalar * vals[src][col];
        }
    }
    public void rowScale(double scalar, int row){
        for(int col=0; col<cols; col++){
            vals[row][col] *= scalar;
        }
    }
    public String toString(){
        StringBuilder sb = new StringBuilder();
        for(int row=0; row<rows; row++){
            for(int col=0; col<cols; col++){
                sb.append(vals[row][col]).append(" ");
            }
            sb.append("\n");
        }
        return sb.toString();
    }
}
