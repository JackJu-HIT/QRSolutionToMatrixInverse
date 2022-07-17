/**
 * @Function:Solution for Matrix Inverse 
 * @Date:2022-07-17 16:00:00
 * @Author:juchunyu
 * @Last modified:juchunyu
 */
#include <stdio.h>
#include <math.h>

#define SIZE  3    


int main(){
    
    double init[SIZE][SIZE]           = {1,2,-1,3,4,-2,5,-4,1};//需要求解的矩阵
    double Q[SIZE][SIZE]              = {0};
    double R[SIZE][SIZE]              = {0};
    double v[SIZE]                    = {0};
    double QT[SIZE][SIZE]             = {0};
    double RInverse[SIZE][SIZE]       = {0};
    double m_in[SIZE][SIZE]           = {0};
    double m_inverse_[SIZE][SIZE]     = {0};
    double m_inverse[SIZE][SIZE]      = {0}; 
    double MatrixResult[SIZE][SIZE]   = {0};                    //矩阵逆求解结果
    int length = SIZE;  
   
   /**进行QR分解**/
    for(int k = 0;k < length;k++){

        /*R(1:k-1,k) = Q(:,1:k-1)’ * A(:,k)*/
        if(k >= 1){
            for(int i = 0;i < k;i++){
                for(int j = 0;j < length;j++){
                    R[i][k] += Q[j][i]*init[j][k];
                }
            }
        }

        /*v = A(:,k) - Q(:,1:k-1) * R(1:k-1,k)*/
        for(int j = 0;j<length;j++){
            if(k < 1){
             v[j] = init[j][0];
            } else {
                    v[j] = init[j][k];
                    for(int g = 0;g < k;g++){
                        v[j] -= R[g][k]*Q[j][g];
                    }
            }
        }

        /*R(k,k) = norm(v)*/
        for(int i = 0;i < length;i++){
            R[k][k] += v[i]*v[i];
        }

        R[k][k] = sqrt(R[k][k]);

        /*Q(:,k) = v / R(k,k)*/
        for(int i = 0;i < length;i++){
            Q[i][k] = v[i]/R[k][k];
        }  
    }
   
    //求解Q的转置
    for(int i = 0;i < length;i++){
       for(int j = 0;j<length;j++){
          QT[j][i] = Q[i][j];
       }
    }

    /**求解R的逆矩阵***/

    //转置
    for(int i = 0;i < length;i++){
        for(int j = 0;j < length;j++){
            m_in[j][i] = R[i][j];
        }
    }
    
   
    for(int j = 0;j < length;j++){
        //求解对角线
        m_inverse_[j][j] = 1/m_in[j][j];
        for(int i = j+1;i < length;i++){
                double temp = 0;

                for(int k = j;k <= i-1;k++){
                    temp += m_in[i][k]*m_inverse_[k][j];
                }

                m_inverse_[i][j] = -temp/m_in[i][i];
        }
    }
  
    //tranpose
    for(int i = 0;i < length;i++){
        for(int j = 0;j < length;j++){
            m_inverse[j][i] = m_inverse_[i][j];  //逆矩阵
        }
    }
 
    /*A = QR => A^(-1) = R^(-1)*Q^T*/

    for(int i = 0;i < length;i++){
        for(int j = 0;j < length;j++){
            for(int k = 0;k < length;k++){
                MatrixResult[i][j] += m_inverse[i][k]*QT[k][j]; 
            }
        }
    }
  
    
    //打印输出结果
    printf("矩阵求逆的结果：\n");
    for(int i = 0;i < SIZE;++i)
    {
        for(int j = 0;j < SIZE;++j)
        {
            printf("%2.4f    ",MatrixResult[i][j]);
        }
        printf("\n");
    }

}