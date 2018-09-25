
class Gauss {
    constructor() {
       this.ready = false;
       this.second = 0.0;
    }

	next(mean, dev) {
		mean = mean == undefined ? 0.0 : mean;
		dev = dev == undefined ? 1.0 : dev;
		
		if (this.ready) {
			this.ready = false;
			return this.second * dev + mean;
		}
		else {
			var u, v, s;
			do {
				u = 2.0 * Math.random() - 1.0;
				v = 2.0 * Math.random() - 1.0;
				s = u * u + v * v;
			} while (s > 1.0 || s == 0.0);
			
			var r = Math.sqrt(-2.0 * Math.log(s) / s);
			this.second = r * u;
			this.ready = true;
			return r * v * dev + mean;
		}
	};
}

function Determinant(m){
    let result=0;
    if (m.length==1){
        return m[0][0];
    }else if(m.length==2){
        return m[0][0] * m[1][1] - m[0][1] * m[1][0];
    }else if(m.length==3){
        return m[0][0] * m[1][1] * m[2][2] + m[0][1] * m[1][2] * m[2][0] + m[0][2] * m[1][0] * m[2][1] - m[2][0] * m[1][1] * m[0][2] - m[1][0] * m[0][1] * m[2][2] - m[0][0] * m[2][1] * m[1][2];
    }else{
        let m1 = new Array(m.size()-1);//массив N-1 x N-1, значения элементов матрицы порядка N-1
        for (let i=0; i<length-1; i++){
            m1[i] = new Array(m.size()-1);
        }
        for (let i=0; i< m.length; i++){
            for (let j=1; j<m.length; j++){
                for (let k=0; k<m.length; k++){
                    if (k<i){
                        m1[j-1][k] = m[j][k];
                    }else if(k>i){
                        m1[j-1][k-1] = m[j][k];
                    }
                }
            }
            result+= Math.pow(-1,i) *m[0][i] * Determinant(m1);
        }
    }
    return result;
}

function InverseMatrix(A) {
    var det = Determinant(A);
    if (det == 0) return false;
    var N = A.length, invA = [];
    for (var i = 0; i < N; i++)
     { invA[i] = [];
       for (var j = 0; j < N; j++)
        { var B = [], sign = ((i+j)%2==0) ? 1 : -1;
          for (var m = 0; m < j; m++)
           { B[m] = [];
             for (var n = 0; n < i; n++)   B[m][n] = A[m][n];
             for (var n = i+1; n < N; n++) B[m][n-1] = A[m][n];
           }
          for (var m = j+1; m < N; m++)
           { B[m-1] = [];
             for (var n = 0; n < i; n++)   B[m-1][n] = A[m][n];
             for (var n = i+1; n < N; n++) B[m-1][n-1] = A[m][n];
           }
          invA[i][j] = sign*Determinant(B)/det;
        }
     }
    return invA;
}

function SumMatrix(A,B) {   
    var m = A.length, n = A[0].length, C = [];
    for (var i = 0; i < m; i++)
     { C[i] = [];
       for (var j = 0; j < n; j++) C[i][j] = A[i][j]+B[i][j];
     }
    return C;
}

function SubMatrix(A,B) {   
    var m = A.length, n = A[0].length, C = [];
    for (var i = 0; i < m; i++)
     { C[i] = [];
       for (var j = 0; j < n; j++) C[i][j] = A[i][j]-B[i][j];
     }
    return C;
}

function MulMatrixNumber(a,A) {   
    var m = A.length, n = A[0].length, B = [];
    for (var i = 0; i < m; i++)
     { B[i] = [];
       for (var j = 0; j < n; j++) B[i][j] = a*A[i][j];
     }
    return B;
}

function MulMatrix(A,B) {
    var rowsA = A.length, colsA = A[0].length,
        rowsB = B.length, colsB = B[0].length,
        C = [];
    if (colsA != rowsB) return false;
    for (var i = 0; i < rowsA; i++) C[i] = [];
    for (var k = 0; k < colsB; k++)
     { for (var i = 0; i < rowsA; i++)
        { var t = 0;
          for (var j = 0; j < rowsB; j++) t += A[i][j]*B[j][k];
          C[i][k] = t;
        }
     }
    return C;
}

function MulMatrixVec(V,B) {
    var rowsA = 1, colsA = V.length,
        rowsB = B.length, colsB = B[0].length,
        C = [];
    if (colsA != rowsB) return false;
    
    for (var k = 0; k < colsB; k++) { 
        C[k] = 0;
        for (var i = 0; i < colsA; i++) { 
            var t = 0;
            C[k] += V[i]*B[i][k];
        }
     }
    return C;
}

function MulVecVec(V1, V2) {
    var C = 0;
    for (let i = 0; i < V1.length; i+=1) {
        C += V1[i]*V2[i];
    }
    return C;
}

function SubVecVec(V1, V2) {
    var C = [];
    for (let i = 0; i < V1.length; i+=1) {
        C.push(V1[i]-V2[i]);
    }
    return C;
}

function TransMatrix(A) {
    var m = A.length, n = A[0].length, AT = [];
    for (var i = 0; i < n; i++)
     { AT[i] = [];
       for (var j = 0; j < m; j++) AT[i][j] = A[j][i];
     }
    return AT;
}


class NGauss {
    constructor() {
        this.n = 0;
        this.rNorm = new Gauss();
        this.X = [];
    }
    setN(_n) {
        this.n = _n;
    }
    setM(_m) {
        this.M = _m;
    }
    setD(_d) {
        this.D = _d;
    }
    setK(_k) {
        this.K = _k;
        let isNorm = document.getElementById('isNorm').checked;
        if (isNorm) {
            for (let i = 0; i < this.n ; i+=1) {
                for (let j = 0; j < this.n ; j+=1) {
                    this.K[i][j] *= Math.sqrt(this.D[i]*this.D[j])
                }
            }
        }
        let det = Determinant(this.K);
        if (det<=0) {
            console.log('PROBLEM',det);
        }
        this.A = new Array(this.n);
        for (let i = 0; i < this.n; i+=1) {
            this.A[i] = new Array(this.n);
        }
        for (let i = 0; i < this.n; i+=1) {
            for (let j = 0; j < this.n; j+=1) {
                if (j <= i) {
                    let sumaijk = 0.0;
                    let sumajk = 0.0;
                    for (let k = 0; k < j; k +=1) {
                        sumaijk += this.A[i][k]*this.A[j][k];
                        sumajk += this.A[j][k]*this.A[j][k];
                    }
                    this.A[i][j] = (this.K[i][j]-sumaijk) / 
                        Math.sqrt(this.K[j][j] - sumajk);
                    //console.log(this.A[i][j], '=', '(',this.K[i][j],'-',sumaijk,') / Math.sqrt(',this.K[j][j], '-', sumajk,')');
                } else {
                    this.A[i][j] = 0;
                }
            }
        }
        //console.log("M",this.M);
        //console.log("D",this.D);
        //console.log("K",this.K);
        //console.log("A",this.A);
    }
    calcParams(inputSize) {
        this.MS = new Array(this.n);
        for (let i = 0; i < this.n; i+=1) {
            this.MS[i] = 0;
        }
        this.DS = new Array(this.n);
        for (let i = 0; i < this.n; i+=1) {
            this.DS[i] = 0;
        }
        this.KS = new Array(this.n);
        for (let i = 0; i < this.n; i+=1) {
            this.KS[i] = new Array(this.n);
        }

        //MA
        for (let i = 0; i < inputSize; i+=1) {
            for (let j = 0; j < this.n; j+=1) {
                this.MS[j] += this.X[i][j];
            }
        }
        for (let j = 0; j < this.n; j+=1) {
            this.MS[j] /= inputSize;
        }
        //DA
        for (let i = 0; i < inputSize; i+=1) {
            for (let j = 0; j < this.n; j+=1) {
                this.DS[j] += ((this.X[i][j]-this.MS[j])*(this.X[i][j]-this.MS[j]));
            }
        }
        for (let j = 0; j < this.n; j+=1) {
            this.DS[j] /= inputSize;
        }
        //KA
        for (let j = 0; j < this.n; j+=1) {
            for (let k = 0; k < this.n; k+=1) {
                this.KS[j][k] = 0;
            }
        }
        for (let i = 0; i < inputSize; i+=1) {
            for (let j = 0; j < this.n; j+=1) {
                for (let k = 0; k < this.n; k+=1) {
                    this.KS[j][k] += ((this.X[i][j]-this.MS[j])*(this.X[i][k]-this.MS[k]));
                }
            }
        }
        for (let j = 0; j < this.n; j+=1) {
            for (let k = 0; k < this.n; k+=1) {
                this.KS[j][k] /= inputSize;
            }
        }
    }
    gen() {
        let X = [];
        let vec = [];
        for (let i = 0; i < this.n; i+=1) {
            vec.push(this.rNorm.next());
        }
        //X = AN+M
        for (let i = 0; i < this.n ; i+=1) {
            let _x = 0;
            for (let j = 0; j < this.n ; j+=1) {
                _x += (this.A[i][j] * vec[j]);
            }
            _x += this.M[i];
            X.push(_x);
        }
        return X;
    }

}

class Mix {
    constructor() {
        this.X = new Array();
        this.Xtrain = new Array();
        this.P = new Array();
        this.Xs = new Array();

        this.C = new Array();
        this.T = new Array();

        this.Perr = 0.0;
    }

    calcParams() {
        this.T = new Array(this.P.length);
        for (let i = 0; i < this.P.length; i+=1) {
            this.T[i] = new Array(this.P.length);
            for (let j = 0; j < this.P.length; j+=1) {
                this.T[i][j] = 0;
            }
        }

        for (let i = 0; i < this.X.length; i+=1) {
            for (let k = 0; k < this.P.length; k+=1) {
                //for (let k = 0; k < this.P.length; k+=1) {

                    let xx = this.X[i].val;
                    let j = this.X[i].pclss;
                    if (j === k) {
                        continue;
                    }
                    let xM1 = SubVecVec(xx, randNormNs[j].M);
                    let xM2 = SubVecVec(xx, randNormNs[k].M); 
                    let r1 = -0.5*MulVecVec(MulMatrixVec(xM1, InverseMatrix(randNormNs[j].K)), xM1) -
                        Math.log(rNormMix.P[j]);    
                    let r2 = -0.5*MulVecVec(MulMatrixVec(xM2, InverseMatrix(randNormNs[k].K)), xM2) -
                        Math.log(rNormMix.P[k]);  
                    
                    this.X[i].pclss = (r1 > r2) ? j : k;

                    //let res12 = MulVecVec(MulMatrixVec(xM1, InverseMatrix(randNormNs[j].K)), xM1) -
                    //    MulVecVec(MulMatrixVec(xM2, InverseMatrix(randNormNs[k].K)), xM2) + 
                    //    Math.log(Determinant(randNormNs[j].K)/Determinant(randNormNs[k].K));
                    //this.X[i].pclss = (res12 <= 2*Math.log(rNormMix.P[j]/rNormMix.P[k])) ? j : k;
                //}
            }
            this.T[this.X[i].clss][this.X[i].pclss] += 1;
        }
        let err = 0.0;
        for (let i = 0; i < this.X.length; i+=1) {
            if (this.X[i].clss !== this.X[i].pclss) {
                err+=1;
            }
        }
        this.Perr = err / this.X.length;
        //console.log(this.X);
        console.log(this.Perr);
        console.log(this.T);
    }

    calcParamsWithMaxEntrop() {
        this.T = new Array(this.P.length);
        for (let i = 0; i < this.P.length; i+=1) {
            this.T[i] = new Array(this.P.length);
            for (let j = 0; j < this.P.length; j+=1) {
                this.T[i][j] = 0;
            }
        }

        for (let i = 0; i < this.X.length; i+=1) {
            for (let k = 0; k < this.P.length; k+=1) {
                //for (let k = 0; k < this.P.length; k+=1) {

                    let xx = this.X[i].val;
                    let j = this.X[i].pclss;
                    if (j === k) {
                        continue;
                    }
                    let xM1 = SubVecVec(xx, randNormNs[j].MS);
                    let xM2 = SubVecVec(xx, randNormNs[k].MS); 
                    let r1 = -0.5*MulVecVec(MulMatrixVec(xM1, InverseMatrix(randNormNs[j].KS)), xM1) -
                        Math.log(rNormMix.P[j]);    
                    let r2 = -0.5*MulVecVec(MulMatrixVec(xM2, InverseMatrix(randNormNs[k].KS)), xM2) -
                        Math.log(rNormMix.P[k]);  
                    
                    this.X[i].pclss = (r1 > r2) ? j : k;

                    //let res12 = MulVecVec(MulMatrixVec(xM1, InverseMatrix(randNormNs[j].K)), xM1) -
                    //    MulVecVec(MulMatrixVec(xM2, InverseMatrix(randNormNs[k].K)), xM2) + 
                    //    Math.log(Determinant(randNormNs[j].K)/Determinant(randNormNs[k].K));
                    //this.X[i].pclss = (res12 <= 2*Math.log(rNormMix.P[j]/rNormMix.P[k])) ? j : k;
                //}
            }
            this.T[this.X[i].clss][this.X[i].pclss] += 1;
        }
        let err = 0.0;
        for (let i = 0; i < this.X.length; i+=1) {
            if (this.X[i].clss !== this.X[i].pclss) {
                err+=1;
            }
        }
        this.Perr = err / this.X.length;
        //console.log(this.X);
        console.log(this.Perr);
        console.log(this.T);
    }


    calcParamsWithParWin() {
        this.T = new Array(this.P.length);
        for (let i = 0; i < this.P.length; i+=1) {
            this.T[i] = new Array(this.P.length);
            for (let j = 0; j < this.P.length; j+=1) {
                this.T[i][j] = 0;
            }
        }

        let ker = (d, u) => {
            return 1/Math.pow(Math.PI*2, d/2) * 
                Math.exp(-0.5*u*u);
        }


        for (let i = 0; i < this.X.length; i+=1) {
            for (let k = 0; k < this.P.length; k+=1) {
                let xx = this.X[i].val;
                let j = this.X[i].pclss;
                if (j === k) {
                    continue;
                }

                let metric = (x, y) => {
                    let res = 0;
                    for (let i = 0; i < x.length; i+=1) {
                        res += (x[i]-y[i])*(x[i]-y[i])
                    }
                    return Math.sqrt(res);
                }

                let to = randNormNs[k].X.length < randNormNs[j].X.length ?
                    randNormNs[k].X.length : randNormNs[j].X.length;

                let r1arr = [];
                for (let l = 0; l < randNormNs[k].X.length; l+=1) {
                    r1arr.push(metric(xx, randNormNs[k].X[l]));
                }
                r1arr.sort((a, b) => {
                    return a - b;
                });
                let r2arr = [];
                for (let l = 0; l < randNormNs[j].X.length; l+=1) {
                    r2arr.push(metric(xx, randNormNs[j].X[l]));
                }
                r2arr.sort((a, b) => {
                    return a - b;
                });

                //console.log(r1arr, r2arr);

                let r1 = 0;
                for (let l = 0; l < Math.sqrt(to); l+=1) {
                    let v = r1arr[l]/Math.sqrt(to);
                    r1 += ker(this.P.length, v);
                }  
                let r2 = 0;
                for (let l = 0; l < Math.sqrt(to); l+=1) {
                    let v = r2arr[l]/Math.sqrt(to);
                    r2 += ker(this.P.length, v);
                } 
                //console.log('!:: ',r1, r2, j, k);
                this.X[i].pclss = (r1 < r2) ? j : k;
            }
            this.T[this.X[i].clss][this.X[i].pclss] += 1;
        }
        let err = 0.0;
        for (let i = 0; i < this.X.length; i+=1) {
            if (this.X[i].clss !== this.X[i].pclss) {
                err+=1;
            }
        }
        this.Perr = err / this.X.length;
        //console.log(this.X);
        console.log(this.Perr);
        console.log(this.T);
    }

    calcParamsNeighbor() {
        this.T = new Array(this.P.length);
        for (let i = 0; i < this.P.length; i+=1) {
            this.T[i] = new Array(this.P.length);
            for (let j = 0; j < this.P.length; j+=1) {
                this.T[i][j] = 0;
            }
        }

        for (let i = 0; i < this.X.length; i+=1) {
            for (let k = 0; k < this.P.length; k+=1) {
                let xx = this.X[i].val;
                let j = this.X[i].pclss;
                if (j === k) {
                    continue;
                }

                let metric = (x, y) => {
                    let res = 0;
                    for (let i = 0; i < x.length; i+=1) {
                        res += (x[i]-y[i])*(x[i]-y[i])
                    }
                    return Math.sqrt(res);
                }


                let r1arr = [];
                for (let l = 0; l < randNormNs[k].X.length; l+=1) {
                    r1arr.push(metric(xx, randNormNs[k].X[l]));
                }
                r1arr.sort((a, b) => {
                    return a - b;
                });
                let r2arr = [];
                for (let l = 0; l < randNormNs[j].X.length; l+=1) {
                    r2arr.push(metric(xx, randNormNs[j].X[l]));
                }
                r2arr.sort((a, b) => {
                    return a - b;
                });


                let r1 = 1/(r1arr[0]*r1arr[0]);
                let r2 = 1/(r2arr[0]*r2arr[0]);
                //console.log('!:: ',r1, r2, j, k);
                this.X[i].pclss = (r1 < r2) ? j : k;
            }
            this.T[this.X[i].clss][this.X[i].pclss] += 1;
        }
        let err = 0.0;
        for (let i = 0; i < this.X.length; i+=1) {
            if (this.X[i].clss !== this.X[i].pclss) {
                err+=1;
            }
        }
        this.Perr = err / this.X.length;
        //console.log(this.X);
        console.log(this.Perr);
        console.log(this.T);
    }

    calcParamsKNeighbor() {
        this.T = new Array(this.P.length);
        for (let i = 0; i < this.P.length; i+=1) {
            this.T[i] = new Array(this.P.length);
            for (let j = 0; j < this.P.length; j+=1) {
                this.T[i][j] = 0;
            }
        }

        for (let i = 0; i < this.X.length; i+=1) {
            for (let k = 0; k < this.P.length; k+=1) {
                let xx = this.X[i].val;
                let j = this.X[i].pclss;
                if (j === k) {
                    continue;
                }

                let metric = (x, y) => {
                    let res = 0;
                    for (let i = 0; i < x.length; i+=1) {
                        res += (x[i]-y[i])*(x[i]-y[i])
                    }
                    return Math.sqrt(res);
                }

                let to = randNormNs[k].X.length < randNormNs[j].X.length ?
                    randNormNs[k].X.length : randNormNs[j].X.length;

                let r1arr = [];
                for (let l = 0; l < randNormNs[k].X.length; l+=1) {
                    r1arr.push(metric(xx, randNormNs[k].X[l]));
                }
                r1arr.sort((a, b) => {
                    return a - b;
                });
                let r2arr = [];
                for (let l = 0; l < randNormNs[j].X.length; l+=1) {
                    r2arr.push(metric(xx, randNormNs[j].X[l]));
                }
                r2arr.sort((a, b) => {
                    return a - b;
                });


                let r1 = 0;
                for (let l = 0; l < Math.sqrt(to); l+=1) {
                    r1 += 1/(r1arr[l]*r1arr[l]);
                }  
                let r2 = 0;
                for (let l = 0; l < Math.sqrt(to); l+=1) {
                    r2 += 1/(r2arr[l]*r2arr[l])
                } 
                //console.log('!:: ',r1, r2, j, k);
                this.X[i].pclss = (r1 < r2) ? j : k;
            }
            this.T[this.X[i].clss][this.X[i].pclss] += 1;
        }
        let err = 0.0;
        for (let i = 0; i < this.X.length; i+=1) {
            if (this.X[i].clss !== this.X[i].pclss) {
                err+=1;
            }
        }
        this.Perr = err / this.X.length;
        //console.log(this.X);
        console.log(this.Perr);
        console.log(this.T);
    }

    calcParamsV() {
        this.T = new Array(this.P.length);
        for (let i = 0; i < this.P.length; i+=1) {
            this.T[i] = new Array(this.P.length);
            for (let j = 0; j < this.P.length; j+=1) {
                this.T[i][j] = 0;
            }
        }

        for (let i = 0; i < this.X.length; i+=1) {
            for (let k = 0; k < this.P.length; k+=1) {
                let xx = this.X[i].val;
                let j = this.X[i].pclss;
                if (j === k) {
                    continue;
                }

                let to = randNormNs[k].X.length < randNormNs[j].X.length ?
                    randNormNs[k].X.length : randNormNs[j].X.length;
                to = Math.sqrt(to);

                let n = 10;
                while (n < to) {
                    n *= 10;
                }
                n /= 10;
                n = 20/to;
                let xxPos = [];
                for (let p = 0; p < xx.length; p+=1) {
                    let f = Math.floor(xx[p]*n)/n;
                    let c = Math.ceil(xx[p]*n)/n;
                    xxPos.push({min: f, max: c});
                }
                //console.log(n, xxPos);

                let testX = (x, pos) => {
                    for (let ii = 0; ii < x.length; ii+=1) {
                        if (x[ii] >= pos[ii].min && x[ii] <= pos[ii].max) {
                            continue;
                        } else {
                            return 0;
                        }
                    }
                    return 1;
                }

                let r1 = 0;
                for (let l = 0; l < to; l+=1) {
                    r1 += testX(randNormNs[k].X[l], xxPos);
                }

                let r2 = 0;
                for (let l = 0; l < to; l+=1) {
                    r2 += testX(randNormNs[j].X[l], xxPos);
                }
                //console.log(r1, r2);
                this.X[i].pclss = (r1 < r2) ? j : k;
            }
            this.T[this.X[i].clss][this.X[i].pclss] += 1;
        }
        let err = 0.0;
        for (let i = 0; i < this.X.length; i+=1) {
            if (this.X[i].clss !== this.X[i].pclss) {
                err+=1;
            }
        }
        this.Perr = err / this.X.length;
        //console.log(this.X);
        console.log(this.Perr);
        console.log(this.T);
    }

    
    calcParamsKNeighborDensity() {
        this.T = new Array(this.P.length);
        for (let i = 0; i < this.P.length; i+=1) {
            this.T[i] = new Array(this.P.length);
            for (let j = 0; j < this.P.length; j+=1) {
                this.T[i][j] = 0;
            }
        }

        for (let i = 0; i < this.X.length; i+=1) {
            for (let k = 0; k < this.P.length; k+=1) {
                let xx = this.X[i].val;
                let j = this.X[i].pclss;
                if (j === k) {
                    continue;
                }

                let metric = (x, y) => {
                    let res = 0;
                    for (let i = 0; i < x.length; i+=1) {
                        res += (x[i]-y[i])*(x[i]-y[i])
                    }
                    return Math.sqrt(res);
                }

                let to = randNormNs[k].X.length < randNormNs[j].X.length ?
                    randNormNs[k].X.length : randNormNs[j].X.length;

                let r1arr = [];
                for (let l = 0; l < randNormNs[k].X.length; l+=1) {
                    r1arr.push({m:metric(xx, randNormNs[k].X[l]), pos:l});
                }
                r1arr.sort((a, b) => {
                    return a.m - b.m;
                });
                let r2arr = [];
                for (let l = 0; l < randNormNs[j].X.length; l+=1) {
                    r2arr.push({m:metric(xx, randNormNs[j].X[l]), pos:l});
                }
                r2arr.sort((a, b) => {
                    return a.m - b.m;
                });

                let getParams = (mix, arr, size) => {
                    let n = randNormNs[0].X[0].length;

                    let MS = new Array(n);
                    for (let i = 0; i < n; i+=1) {
                        MS[i] = 0;
                    }
                    let KS = new Array();
                    for (let i = 0; i < n; i+=1) {
                        KS[i] = new Array(n);
                    }

                    //MS
                    for (let i = 0; i < size; i+=1) {
                        for (let j = 0; j < n; j+=1) {
                            //console.log(randNormNs[mix].X[arr[i].pos][j], '=', arr[i]);
                            MS[j] += randNormNs[mix].X[arr[i].pos][j];
                        }
                    }
                    for (let j = 0; j < n; j+=1) {
                        MS[j] /= size;
                    }

                    for (let j = 0; j < n; j+=1) {
                        for (let k = 0; k < n; k+=1) {
                            KS[j][k] = 0;
                        }
                    }
                    for (let i = 0; i < size; i+=1) {
                        for (let j = 0; j < n; j+=1) {
                            for (let k = 0; k < n; k+=1) {
                                KS[j][k] += (randNormNs[mix].X[arr[i].pos][j] - MS[j]) *
                                    (randNormNs[mix].X[arr[i].pos][k] - MS[k]);
                            }
                        }
                    }
                    for (let j = 0; j < n; j+=1) {
                        for (let k = 0; k < n; k+=1) {
                            KS[j][k] /= size;
                        }
                    }

                    return {K:KS, M:MS};
                }

                let KM1 = getParams(k, r1arr, Math.sqrt(to));

                let KM2 = getParams(j, r2arr, Math.sqrt(to))

                //console.log(KM1, KM2);

                let xM1 = SubVecVec(xx, KM1.M);
                let xM2 = SubVecVec(xx, KM2.M); 
                let r1 = -0.5*MulVecVec(MulMatrixVec(xM1, InverseMatrix(KM1.K)), xM1) -
                    Math.log(0.5);    
                let r2 = -0.5*MulVecVec(MulMatrixVec(xM2, InverseMatrix(KM2.K)), xM2) -
                    Math.log(0.5);  
                
                this.X[i].pclss = (r1 < r2) ? j : k;

            }
            this.T[this.X[i].clss][this.X[i].pclss] += 1;
        }
        let err = 0.0;
        for (let i = 0; i < this.X.length; i+=1) {
            if (this.X[i].clss !== this.X[i].pclss) {
                err+=1;
            }
        }
        this.Perr = err / this.X.length;
        //console.log(this.X);
        console.log(this.Perr);
        console.log(this.T);
    }
}




var rNormN = new NGauss();
var rNormMix = new Mix();

function DrawPointsXs(cl, dim1, dim2, name, isTrain) {//,color, name
    let data = [];
    trace = {
        x: [],
        y: [],
        type: 'scatter',
        mode: 'markers',
        marker: {
            //color: color,
            size: 5,
            symbol: 'circle',
        },
        name: name + " class " + cl,
    };
    
    tracerr = {
        x: [],
        y: [],
        type: 'scatter',
        mode: 'markers',
        marker: {
            //color: color,
            size: 5,
            symbol: 'circle-open',
        },
        name: name + " class error " + cl,
    };

    if (isTrain) {
        for (let i = 0; i < rNormMix.Xtrain.length; i+=1) {
            if (cl == rNormMix.Xtrain[i].clss) {
                let X = rNormMix.Xtrain[i].val[dim1];
                let Y = rNormMix.Xtrain[i].val[dim2];
                trace.x.push(X);
                trace.y.push(Y);
            }
        }
    } else {
        for (let i = 0; i < rNormMix.X.length; i+=1) {
            if (cl == rNormMix.X[i].pclss) {
                if (rNormMix.X[i].clss === rNormMix.X[i].pclss) {
                    let X = rNormMix.X[i].val[dim1];
                    let Y = rNormMix.X[i].val[dim2];
                    trace.x.push(X);
                    trace.y.push(Y);
                } else {
                    let X = rNormMix.X[i].val[dim1];
                    let Y = rNormMix.X[i].val[dim2];
                    tracerr.x.push(X);
                    tracerr.y.push(Y);
                }
            }
        }
    }
    data.push(trace);
    data.push(tracerr);
    return data;
}

var randNormNs = [];
var D = [];


function SeporateFunc() {
    //D = [];
    //let k = 1;
    //let l = 0;

    let data = [];
    
    for (let l = 0; l < rNormMix.P.length; l+=1)
    for (let k = 0; k < rNormMix.P.length; k+=1) {
        if (k <= l) {
            continue;
        }

        let trace = {
            x: [],
            y: [],
            name: 'Seporate func' + l + "_" + k,
            type: 'scatter',
            //color: 'rgb(255, 182, 193)',
            //name: name,
            mode: 'markers',
            marker: {
                //color: color,
                size: 4,
                symbol: 'circle',
            },
        };

        //let xs = [];

        for (let i = -100; i < 100; i+=0.2) 
        for (let j = -100; j < 100; j+=0.2){
            let xx = [i, j];
            
            let Bkl = SubMatrix(InverseMatrix(randNormNs[k].K), InverseMatrix(randNormNs[l].K));
            let xBkl = MulMatrixVec(xx, Bkl);
            let xBklx = MulVecVec(xBkl, xx);
            let MB = 2*MulVecVec(SubVecVec(MulMatrixVec(randNormNs[l].M, InverseMatrix(randNormNs[l].K)),
                MulMatrixVec(randNormNs[k].M, InverseMatrix(randNormNs[k].K))), xx);
            
            let LN = Math.log(Determinant(randNormNs[l].K)/Determinant(randNormNs[k].K)) +
                2*Math.log(rNormMix.P[l]/rNormMix.P[k]) -
                MulVecVec(MulMatrixVec(randNormNs[l].M, InverseMatrix(randNormNs[l].K)), randNormNs[l].M) +
                MulVecVec(MulMatrixVec(randNormNs[k].M, InverseMatrix(randNormNs[k].K)), randNormNs[k].M);
            if (xBklx + MB + LN <= 0.2 && xBklx + MB + LN >= -0.2) {
                
                //xs.push({x:xx[0], y:xx[1]});
                trace.x.push(xx[0]);
                trace.y.push(xx[1]);
            }
            
            //let xM1 = SubVecVec(xx, randNormNs[l].M);
            //let xM2 = SubVecVec(xx, randNormNs[k].M); 
            //let r1 = -0.5*MulVecVec(MulMatrixVec(xM1, InverseMatrix(randNormNs[l].K)), xM1) -
            //    Math.log(rNormMix.P[l]);    
            //let r2 = -0.5*MulVecVec(MulMatrixVec(xM2, InverseMatrix(randNormNs[k].K)), xM2) -
            //    Math.log(rNormMix.P[k]);  
            //
            //if (r1-r2 <= 1 && r1-r2 >= -1) {
            //    trace.x.push(xx[0]);
            //    trace.y.push(xx[1]);
            //    //trace.y.push(xBklx + MB + LN);
            //}
        }  
        data.push(trace);
    }
    return data;
}

function randNormN() {
    randNormNs = [];
    let inputMix = parseInt(document.getElementById('inputMix').value);
    let inputSize = parseInt(document.getElementById('inputSize').value);
    let inputSizeTrain = parseInt(document.getElementById('inputSizeTrain').value);
    let inputN = parseInt(document.getElementById('inputN').value);
    rNormMix.Xs = [];
    for (let k = 0; k < inputMix; k+=1) {
        let _rNormN = new NGauss();
        let D = new Array(inputN);
        let M = new Array(inputN);
        let K = new Array(inputN);
        for (let i = 0; i < inputN; i+=1) {
            K[i] = new Array(inputN);
        }
        
        for (let i = 0; i < inputN; i+=1) {
            for (let j = 0; j < inputN; j+=1) {
                K[i][j] = parseFloat(document.getElementById('inK'+k+(i*inputN+j)).value);
            }
            D[i] = parseFloat(document.getElementById('inD'+k+i).value);
            M[i] = parseFloat(document.getElementById('inM'+k+i).value);
        }
        _rNormN.setN(inputN);
        _rNormN.setM(M);
        _rNormN.setD(D);
        _rNormN.setK(K);
        
        _rNormN.X = [];

        randNormNs.push(_rNormN);
    }

    rNormMix.P = [];
    for (let i = 0; i < inputMix; i+=1) {
        rNormMix.P.push(parseFloat(document.getElementById('inP'+i).value)); 
    }

    getClss = () => {
        let r = Math.random();

        let p = rNormMix.P[0];
        let clss = 0;
        for (let j=1;j<rNormMix.P.length; j+=1) {
            if (p < r) {
                p += rNormMix.P[j];
                clss += 1; 
            }
        }
        return clss;
    }

    rNormMix.X = []
    for (let i = 0; i < inputSize; i+=1) {
        let clss = getClss();
        rNormMix.X.push({
            val: randNormNs[clss].gen(),
            clss: clss,
            pclss: 0,
        });
    }

    rNormMix.Xtrain = []
    for (let i = 0; i < randNormNs.length; i+=1) {
        randNormNs[i].X = [];
    }
    for (let i = 0; i < inputSizeTrain; i+=1) {
        let clss = getClss();

        let val = randNormNs[clss].gen();

        randNormNs[clss].X.push(val);
        //randNormNs[clss].setN(randNormNs[clss].X.length);

        rNormMix.Xtrain.push({
            val: val,
            clss: clss,
            pclss: 0,
        });
    }

    for (let i = 0; i < randNormNs.length; i+=1) {
        randNormNs[i].calcParams(randNormNs[i].X.length);
    }

    console.log('Mix ',rNormMix);
    console.log('Ns ',randNormNs);

    let time = performance.now();

    let listval = document.getElementById('blist').value;
    switch(listval) {
        case 'Идеальный к-р':
            rNormMix.calcParams();
        break;
        case 'Максимального правдоподобия':
            rNormMix.calcParamsWithMaxEntrop();
        break;
        case 'Парзеновские окона':
            rNormMix.calcParamsWithParWin();
        break;
        case 'Ближайший сосед':
            rNormMix.calcParamsNeighbor();
        break;
        case 'K ближайших соседей':
            rNormMix.calcParamsKNeighbor();    
        break;
        case 'Оценка плотности через K ближайших соседей':
            rNormMix.calcParamsKNeighborDensity();    
        break;
        case 'Через ячейки фиксированного объёма':
            rNormMix.calcParamsV();    
        break;
        default:
        break;
    }
    time = performance.now() - time;
    printParams(time);
    ////////DRAW
    Draw();
}

function Draw() {
    let inputMix = parseInt(document.getElementById('inputMix').value);
    let inputPoints = document.getElementById('printPoints').checked;
    let inputPointsTrain = document.getElementById('printPointsTrain').checked;
    let inputFunc = document.getElementById('printFunc').checked;
    let data = [];

    //console.log(inputPointsTrain);

    for (let i = 0; i < inputMix; i+=1) {
        if (inputPoints) {
            let __pp = DrawPointsXs(
                i, dim1, dim2, '', false 
            );
            data.push(__pp[0]);
            data.push(__pp[1]);
        }
        if (inputPointsTrain) {
            data.push(DrawPointsXs(
                i, dim1, dim2, 'train', true 
            )[0]);
        }
    }
    if (inputFunc) {
        let listval = document.getElementById('blist').value;
        switch(listval) {
            case 'Идеальный к-р':
            
            break;
            case 'Максимального правдоподобия':
            
            break;
            case 'Парзеновские окона':
            
            break;
            default:
            break;  
        }
        let _d = SeporateFunc();
        for (let o = 0; o < _d.length; o+=1) {
            data.push(_d[o]);
        }
    }
        
    let layout = {}
    layout = {
        title: 'Visualization',
        width: 700,
        height: 700,

        yaxis:{
            autotick: true,
            scaleanchor: "x",
        },
        xaxis:{
            autotick: true,
        },
    };
    
    Plotly.newPlot('graph', data, layout);
}