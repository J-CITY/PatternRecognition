
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

        this.NewClasses = new Array();
        //for T alg
        this.T = 0;
        //for Max dist
        this.Alpha = 0;
        //for KMid
        this.ClassSize = 0;
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

    calcParamsTalg() {
        this.NewClasses = new Array();
        this.T = parseFloat(document.getElementById('T').value);

        let metric = (x, y) => {
            let res = 0;
            for (let i = 0; i < x.length; i+=1) {
                res += (x[i]-y[i])*(x[i]-y[i])
            }
            return Math.sqrt(res);
        }

        this.NewClasses.push({X: [], M: []});
        this.NewClasses[0].X.push({
            val:rNormMix.Xtrain[0].val,
            clss:rNormMix.Xtrain[0].clss,
        });
        this.NewClasses[0].M = rNormMix.Xtrain[0].val;

        for (let i = 1; i < rNormMix.Xtrain.length; i+=1) {
            let _x = rNormMix.Xtrain[i].val;
            
            let minMetric = 1000000;
            let clss = 0;
            for (let j = 0; j < this.NewClasses.length; j+=1) {
                let metr = metric(_x, this.NewClasses[j].M);
                if (metr < this.T && metr < minMetric) {
                    minMetric = metr;
                    clss = j;
                }
            }

            if (minMetric === 1000000) {
                this.NewClasses.push({X: [], M: []});
                this.NewClasses[this.NewClasses.length-1].X.push({
                    val: rNormMix.Xtrain[i].val,
                    clss: rNormMix.Xtrain[i].clss,
                });
                this.NewClasses[this.NewClasses.length-1].M = rNormMix.Xtrain[i].val;
            } else {
                this.NewClasses[clss].X.push({
                    val: rNormMix.Xtrain[i].val,
                    clss: rNormMix.Xtrain[i].clss,
                });
                //New M
                let n = _x.length;
                let MS = new Array(n);
                for (let k = 0; k < n; k+=1) {
                    MS[k] = 0;
                }
                for (let k = 0; k < this.NewClasses[clss].X.length; k+=1) {
                    for (let kk = 0; kk < n; kk+=1) {
                        MS[kk] += this.NewClasses[clss].X[k].val[kk];
                    }
                }
                for (let k = 0; k < n; k+=1) {
                    MS[k] /= this.NewClasses[clss].X.length;
                }
                this.NewClasses[clss].M = MS;
            }
        }
    }

    calcParamsMaxDist() {
        this.Alpha = parseFloat(document.getElementById('Alpha').value);

        this.NewClasses = new Array();

        let metric = (x, y) => {
            let res = 0;
            for (let i = 0; i < x.length; i+=1) {
                res += (x[i]-y[i])*(x[i]-y[i])
            }
            return Math.sqrt(res);
        }
        //z1
        this.NewClasses.push({X: [], M: []});
        //this.NewClasses[0].X.push(rNormMix.Xtrain[0]);
        this.NewClasses[0].M = rNormMix.Xtrain[0].val;

        //z2
        let maxMetric = 0;
        let maxMetricId = 0;
        for (let i = 1; i < rNormMix.Xtrain.length; i+=1) {
            let metr = metric(rNormMix.Xtrain[0].val, rNormMix.Xtrain[i].val);
            if (metr > maxMetric) {
                maxMetric = metr;
                maxMetricId = i;
            }
        }
        this.NewClasses.push({X: [], M: []});
        //this.NewClasses[1].X.push(rNormMix.Xtrain[maxMetricId]);
        this.NewClasses[1].M = rNormMix.Xtrain[maxMetricId].val;

        //zi
        let test = true;
        while (test) {
            let maxMetric = 0;
            let maxMetricId = 0;
            for (let i = 1; i < rNormMix.Xtrain.length; i+=1) {
                let minMetric = 1000000;
                let minMetricId = 0;
                for (let j = 0; j < this.NewClasses.length; j+=1) {
                    let metr = metric(this.NewClasses[j].M, rNormMix.Xtrain[i].val);
                    if (metr < minMetric) {
                        minMetric = metr;
                        minMetricId = i;
                    }
                }
                if (maxMetric < minMetric) {
                    maxMetric = minMetric;
                    maxMetricId = minMetricId;
                }
            }
            for (let j = 0; j < this.NewClasses.length; j+=1) {
                let metr = metric(this.NewClasses[j].M, rNormMix.Xtrain[maxMetricId].val);
                if (metr*this.Alpha >= maxMetric) {
                    test = false;
                    continue;
                }
            }

            this.NewClasses.push({X: [], M: []});
            //this.NewClasses[this.NewClasses.length-1].X.push(rNormMix.Xtrain[maxMetricId]);
            this.NewClasses[this.NewClasses.length-1].M = rNormMix.Xtrain[maxMetricId].val;
        }
            
        for (let i = 0; i < rNormMix.Xtrain.length; i+=1) {
            let _x = rNormMix.Xtrain[i].val;

            let minMetric = 1000000;
            let clss = 0;
            for (let j = 0; j < this.NewClasses.length; j+=1) {
                let metr = metric(_x, this.NewClasses[j].M);
                //console.log(metr, i, _x);
                if (metr < minMetric) {
                    minMetric = metr;
                    clss = j;
                }
            }

            this.NewClasses[clss].X.push({
                val: rNormMix.Xtrain[i].val,
                clss: rNormMix.Xtrain[i].clss,
            });
            //New M
            let n = _x.length;
            let MS = new Array(n);
            for (let k = 0; k < n; k+=1) {
                MS[k] = 0;
            }
            for (let k = 0; k < this.NewClasses[clss].X.length; k+=1) {
                for (let kk = 0; kk < n; kk+=1) {
                    MS[kk] += this.NewClasses[clss].X[k].val[kk];
                }
            }
            for (let k = 0; k < n; k+=1) {
                MS[k] /= this.NewClasses[clss].X.length;
            }
            this.NewClasses[clss].M = MS;
        }
    }

    calcParamsKMid() {
        this.ClassSize = parseFloat(document.getElementById('ClassCount').value);
        this.eps = parseInt(document.getElementById('eps').value);


        this.NewClasses = new Array();

        let metric = (x, y) => {
            let res = 0;
            for (let i = 0; i < x.length; i+=1) {
                res += (x[i]-y[i])*(x[i]-y[i])
            }
            return Math.sqrt(res);
        }
        //z1...zi
        for (let i = 0; i < this.ClassSize; i+=1) {
            this.NewClasses.push({
                X: [], 
                M: rNormMix.Xtrain[i].val,
            });
        }

        let test = true;
        while (test) {
            for (let j = 0; j < this.NewClasses.length; j+=1) {
                this.NewClasses[j].X = [];
            }
            for (let i = 0; i < rNormMix.Xtrain.length; i+=1) {
                let minMetric = 1000000;
                let minMetricId = 0;
                let clss = 0;
                for (let j = 0; j < this.NewClasses.length; j+=1) {
                    let metr = metric(this.NewClasses[j].M, rNormMix.Xtrain[i].val);
                    if (metr < minMetric) {
                        minMetric = metr;
                        minMetricId = i;
                        clss = j;
                    }
                }
                this.NewClasses[clss].X.push({
                    val: rNormMix.Xtrain[minMetricId].val,
                    clss: rNormMix.Xtrain[minMetricId].clss,
                });
            }


            let epsTest = (x, y, eps) => {
                for (let i = 0; i < x.length; i+=1) {
                    if (Math.abs(x[i] - y[i]) > eps) {
                        return false;
                    }
                }
                return true;
            }

            let _t = true;
            for (let i = 0; i < this.NewClasses.length; i+=1) {
                //New M
                let n = this.NewClasses[0].X[0].val.length;
                let MS = new Array(n);
                for (let k = 0; k < n; k+=1) {
                    MS[k] = 0;
                }
                for (let k = 0; k < this.NewClasses[i].X.length; k+=1) {
                    for (let kk = 0; kk < n; kk+=1) {
                        MS[kk] += this.NewClasses[i].X[k].val[kk];
                    }
                }
                for (let k = 0; k < n; k+=1) {
                    MS[k] /= this.NewClasses[i].X.length;
                }

                _t = _t && epsTest(this.NewClasses[i].M, MS, this.eps);
                this.NewClasses[i].M = MS;
            }
            if (_t) {
                test = false;
            }
        }   
    }

    _calcParams() {
        this.T = new Array(this.NewClasses.length);
        for (let i = 0; i < this.T.length; i+=1) {
            this.T[i] = new Array(this.P.length);
        }
        for (let i = 0; i < this.T.length; i+=1) {
            for (let j = 0; j < this.T[i].length; j+=1) {
                this.T[i][j] = 0;
            }
        }
        for (let i = 0; i < this.NewClasses.length; i+=1) {
            for (let j = 0; j < this.NewClasses[i].X.length; j+=1) {
                this.T[i][this.NewClasses[i].X[j].clss]+=1;
            }
        }
    }

    
    calcParamsIsodata() {
        this.ClassSize = parseInt(document.getElementById('ClassCount').value);
        this.Nc = parseInt(document.getElementById('Nc').value);
        this.eps = parseFloat(document.getElementById('eps').value);
        this.On = parseFloat(document.getElementById('On').value);
        this.Os = parseFloat(document.getElementById('Os').value);
        this.Oc = parseFloat(document.getElementById('Oc').value);
        this.L = parseFloat(document.getElementById('L').value);
        this.I = parseFloat(document.getElementById('I').value);
        this.k = parseFloat(document.getElementById('k').value);


        this.NewClasses = new Array();

        let metric = (x, y) => {
            let res = 0;
            for (let i = 0; i < x.length; i+=1) {
                res += (x[i]-y[i])*(x[i]-y[i])
            }
            return Math.sqrt(res);
        }

        let n11 = () => {
            //11
            let DD = [];
            for (let i = 0; i < this.NewClasses.length; i+=1) {
                for (let j = 0; j < this.NewClasses.length; j+=1) {
                    if (j > i) {
                        let d = metric(this.NewClasses[i].M, this.NewClasses[j].M);
                        if (d < this.Oc) {
                            DD.push({
                                d: d,
                                i: i,
                                j: j,
                            })
                        }
                    }
                }
            }
            DD.sort((a, b) => {
                return a.d - b.d;
            });
            let _set = new Set();
            let _L=0;
            for (let i = 0; i < DD.length; i+=1) {
                if (!_set.has(DD[i].i) && !_set.has(DD[i].j) && _L < this.L) {
                    _L+=1;
                    _set.add(DD[i].i);
                    _set.add(DD[i].j);

                    let newM = [];
                    for (let j = 0; j < this.NewClasses[DD[i].i].M.length; j+=1) {
                        newM.push((1/(this.NewClasses[DD[i].i].X.length + this.NewClasses[DD[i].j].X.length))*
                            (this.NewClasses[DD[i].i].X.length*this.NewClasses[DD[i].i].M[j] +
                            this.NewClasses[DD[i].j].X.length*this.NewClasses[DD[i].j].M[j]));
                    }

                    this.NewClasses.push({
                        X: [], 
                        M: newM,
                        Dm: 0,
                        Sigma: [],
                        SigmaMax: 0,
                        SigmaMaxId: 0,
                    });
                    this.Nc -= 1;
                }
            }
            let nc = [];
            for (let i = 0; i < this.NewClasses.length; i+=1) {
                if (_set.has(i)) {
                    continue;
                } else {
                    nc.push(this.NewClasses[i]);
                }
            }
            this.NewClasses = nc;
        }


        //1
        for (let i = 0; i < this.Nc; i+=1) {
            this.NewClasses.push({
                X: [], 
                M: rNormMix.Xtrain[i].val,
                Dm: 0,
                Sigma: [],
                SigmaMax: 0,
                SigmaMaxId: 0,
            });
        }

        let test = true;
        let iter = 0;
        while (test) {
            //2
            for (let j = 0; j < this.NewClasses.length; j+=1) {
                this.NewClasses[j].X = [];
            }
            for (let i = 0; i < rNormMix.Xtrain.length; i+=1) {
                let minMetric = 1000000;
                let minMetricId = 0;
                let clss = 0;
                for (let j = 0; j < this.NewClasses.length; j+=1) {
                    let metr = metric(this.NewClasses[j].M, rNormMix.Xtrain[i].val);
                    if (metr < minMetric) {
                        minMetric = metr;
                        minMetricId = i;
                        clss = j;
                    }
                }
                this.NewClasses[clss].X.push({
                    val: rNormMix.Xtrain[minMetricId].val,
                    clss: rNormMix.Xtrain[minMetricId].clss,
                });
            }
            //3
            for (let j = 0; j < this.NewClasses.length; j+=1) {
                if (this.NewClasses[j].X.length < this.On) {
                    this.NewClasses.splice(this.NewClasses.indexOf(j), 1);
                    this.Nc -= 1;
                }
            }
            //4
            for (let i = 0; i < this.NewClasses.length; i+=1) {
                //New M
                let n = this.NewClasses[0].X[0].val.length;
                let MS = new Array(n);
                for (let k = 0; k < n; k+=1) {
                    MS[k] = 0;
                }
                for (let k = 0; k < this.NewClasses[i].X.length; k+=1) {
                    for (let kk = 0; kk < n; kk+=1) {
                        MS[kk] += this.NewClasses[i].X[k].val[kk];
                    }
                }
                for (let k = 0; k < n; k+=1) {
                    MS[k] /= this.NewClasses[i].X.length;
                }
                //console.log(MS, this.NewClasses[i]);
                this.NewClasses[i].M = MS;
            }
            //5
            for (let i = 0; i < this.NewClasses.length; i+=1) {
                //New Ds
                let n = this.NewClasses[0].X[0].val.length;
                let DS = 0;
                let d = 0;
                for (let k = 0; k < this.NewClasses[i].X.length; k+=1) {
                    d += metric(this.NewClasses[i].X[k].val, this.NewClasses[i].M);
                }
                d /= this.NewClasses[i].X.length;
                this.NewClasses[i].Ds = d;
            }
            //6
            let DsMid = 0
            for (let i = 0; i < this.NewClasses.length; i+=1) {
                //New Ds
                DsMid += this.NewClasses[i].Ds*this.NewClasses[i].X.length;
            }
            DsMid /= rNormMix.Xtrain.length;

            //7
            iter += 1;

            if (iter == this.I) {
                this.Oc = 0;
                n11();
                break;
            } else
            if (this.Nc >= this.ClassSize*2 || iter % 2 === 0) {
                n11();
                break;
            } else 
            if (this.Nc <= this.ClassSize / 2) {
                //8
                for (let k = 0; k < this.NewClasses.length; k+=1) {
                    let n = this.NewClasses[0].X[0].val.length;
                    let DS = new Array(n);
                    for (let i = 0; i < n; i+=1) {
                        DS[i] = 0;
                    }
                    for (let i = 0; i < this.NewClasses[k].X.length; i+=1) {
                        for (let j = 0; j < n; j+=1) {
                            DS[j] += ((this.NewClasses[k].X[i].val[j]-this.NewClasses[k].M[j])*
                                (this.NewClasses[k].X[i].val[j]-this.NewClasses[k].M[j]));
                        }
                    }  
                    for (let j = 0; j < n; j+=1) {
                        DS[j] /= this.NewClasses[k].X.length;
                        DS[j] = Math.sqrt(DS[j]);
                    }
                    this.NewClasses[k].Sigma = DS;
                }
                //9
                for (let k = 0; k < this.NewClasses.length; k+=1) {
                    let maxSigma = 0;
                    let maxSigmaId = 0;
                    for (let i = 0; i < this.NewClasses[k].Sigma.length; i+=1) {
                        if (maxSigma < this.NewClasses[k].Sigma[i]) {
                            maxSigma = this.NewClasses[k].Sigma[i];
                            maxSigmaId = i;
                        }
                    }  
                    this.NewClasses[k].SigmaMax = maxSigma;
                    this.NewClasses[k].SigmaMaxId = maxSigmaId;
                }
                //10
                for (let k = 0; k < this.NewClasses.length; k+=1) {
                    let jmax = true;
                    if (this.NewClasses[k].SigmaMax <= this.Os) {
                        jmax = false;
                    }
                    let a = true;
                    if (this.NewClasses[k].Ds > DsMid && this.NewClasses[k].X.length > 2*(this.On+1)) {
                        //
                    } else {
                        a = false;
                    }
                    if (jmax && (this.Nc <= this.ClassSize/2 || a)) {
                        let newM1 = this.NewClasses[k].M;
                        newM1[this.NewClasses[k].SigmaMaxId] += this.k * this.NewClasses[k].SigmaMax;
                        this.NewClasses.push({
                            X: [], 
                            M: newM1,
                            Dm: 0,
                            Sigma: [],
                            SigmaMax: 0,
                            SigmaMaxId: 0,
                        });

                        let newM2 = this.NewClasses[k].M;
                        newM2[this.NewClasses[k].SigmaMaxId] -= this.k * this.NewClasses[k].SigmaMax;
                        this.NewClasses.push({
                            X: [], 
                            M: newM2,
                            Dm: 0,
                            Sigma: [],
                            SigmaMax: 0,
                            SigmaMaxId: 0,
                        });
                        this.NewClasses.splice(this.NewClasses.indexOf(k), 1);
                        this.Nc += 1;
                    }
                }
                //n11();
            }
        }
        for (let j = 0; j < this.NewClasses.length; j+=1) {
            this.NewClasses[j].X = [];
        }
        for (let i = 0; i < rNormMix.Xtrain.length; i+=1) {
            let minMetric = 1000000;
            let minMetricId = 0;
            let clss = 0;
            for (let j = 0; j < this.NewClasses.length; j+=1) {
                let metr = metric(this.NewClasses[j].M, rNormMix.Xtrain[i].val);
                if (metr < minMetric) {
                    minMetric = metr;
                    minMetricId = i;
                    clss = j;
                }
            }
            this.NewClasses[clss].X.push({
                val: rNormMix.Xtrain[minMetricId].val,
                clss: rNormMix.Xtrain[minMetricId].clss,
            });
        }

    }
    
}




var rNormN = new NGauss();
var rNormMix = new Mix();
/*
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
*/
function DrawPointsXsTrain(cl, dim1, dim2, name, isTrainTruth) {//,color, name
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
        name: name + "Train class " + cl,
    };


    if (isTrainTruth) {
        for (let i = 0; i < rNormMix.Xtrain.length; i+=1) {
            if (cl == rNormMix.Xtrain[i].clss) {
                let X = rNormMix.Xtrain[i].val[dim1];
                let Y = rNormMix.Xtrain[i].val[dim2];
                trace.x.push(X);
                trace.y.push(Y);
            }
        }
    } else {
        for (let i = 0; i < rNormMix.NewClasses[cl].X.length; i+=1) {
            let X = rNormMix.NewClasses[cl].X[i].val[dim1];
            let Y = rNormMix.NewClasses[cl].X[i].val[dim2];
            trace.x.push(X);
            trace.y.push(Y);
                
        }
    }
    data.push(trace);
    return data;
}

var randNormNs = [];
var D = [];

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

    /*
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
    */
    let listval = document.getElementById('blistauto').value;
    switch(listval) {
        case 'Т алгоритм':
            rNormMix.calcParamsTalg();
        break;
        case 'алгоритм максиминного расстояния':
            rNormMix.calcParamsMaxDist();
        break;
        case 'алгоритм К внутригрупповых средних':
            rNormMix.calcParamsKMid();
        break;
        case 'isodata':
            rNormMix.calcParamsIsodata();
        break
        default:
        break;
    }
    rNormMix._calcParams();
    console.log("New classes", rNormMix.NewClasses);
    printParams();
    ////////DRAW
    Draw();
}

function reDraw() {
    let listval = document.getElementById('blistauto').value;
    switch(listval) {
        case 'Т алгоритм':
            rNormMix.calcParamsTalg();
        break;
        case 'алгоритм максиминного расстояния':
            rNormMix.calcParamsMaxDist();
        break;
        case 'алгоритм К внутригрупповых средних':
            rNormMix.calcParamsKMid();
        break;
        case 'isodata':
            rNormMix.calcParamsIsodata();
        break
        default:
        break;
    }
    rNormMix._calcParams();
    console.log("New classes", rNormMix.NewClasses);
    printParams();
    Draw();
}

function Draw() {
    let inputMix = parseInt(document.getElementById('inputMix').value);
    let inputPointsTrainTruth = document.getElementById('printPointsTrain').checked;
    let inputPointsTrain = document.getElementById('printPointsTrainAuto').checked;
    let data = [];

    //console.log(inputPointsTrain);
    if (inputPointsTrainTruth) {
        for (let i = 0; i < inputMix; i+=1) {
            data.push(DrawPointsXsTrain(
                i, dim1, dim2, ' truth', true 
            )[0]);
        }
    }
    if (inputPointsTrain) {
        for (let i = 0; i < rNormMix.NewClasses.length; i+=1) {
            data.push(DrawPointsXsTrain(
                i, dim1, dim2, '', false 
            )[0]);
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