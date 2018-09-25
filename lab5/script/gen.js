
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
        //console.log("MS",this.MS);
        //console.log("DS",this.DS);
        //console.log("KS",this.KS);

        //console.log('!!!!!!!!!!!!!!!\n');
        //console.log(this.M[0],this.M[1],this.MS[0],this.MS[1],this.M[0]-this.MS[0],this.M[1]-this.MS[1],
        //            this.D[0],this.D[1],this.DS[0],this.DS[1],this.D[0]-this.DS[0],this.D[1]-this.DS[1]
        //);

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


function DrawIsoline(sigma1, sigma2, r12, a1, a2, count, color, mode, name) {
    //console.log(sigma1, sigma2, r12, a1, a2);
    let is3D = document.getElementById('print3D').checked;
    let lines_count = count;
    let trace = {
        x: [],
        y: [],
        z: [],
        type: 'mesh3d',
        mode: mode,
        color: color,
        opacity: 0.5,
        marker: {
            color: color,
            size: 2,
            opacity: 0.5,
        },
        name: name,
    };
    if (!is3D) {
        trace = {
            x: [],
            y: [],
            //type: 'scatter',
            type: 'density',
            
            //color: color,
            mode: mode,
            marker: {
                color: color,
                size: 2,
                opacity: 0.4
            },
            name: name,
        };
    }
    let data = [];

    let chisquareMax = 2.4477;
    let chisquareMin = 0.010;
    let step = (chisquareMax- chisquareMin) / lines_count;

    let chisquare = chisquareMax;
    let discr = Math.sqrt((sigma1*sigma1+sigma2*sigma2)*(sigma1*sigma1+sigma2*sigma2) -
        4*(sigma1*sigma1*sigma2*sigma2-r12*r12));
    let l1 = ((sigma1*sigma1+sigma2*sigma2) + discr) / 2;
    let l2 = ((sigma1*sigma1+sigma2*sigma2) - discr) / 2;

    let angle = -Math.atan2(2*r12*sigma1*sigma2, (sigma1*sigma1-sigma2*sigma2))/2;
    //angle += Math.PI
    //angle /= 2
    
    
    //if (angle < 0) {
    //    angle+= Math.PI*2;
    //}
    //angle += 180
    console.log(angle, "::ANGEL")
    //angle = Math.PI/4;
    //console.log(angle, "::ANGEL")

    while (chisquare > chisquareMin) {
        let halfmajoraxissize = chisquare * Math.sqrt(l1);
        let halfminoraxissize = chisquare * Math.sqrt(l2);


        let alpha = 0;
        while (alpha < 2 * Math.PI) {
            let _X = halfmajoraxissize * Math.cos(alpha);
            let _Y = halfminoraxissize * Math.sin(alpha);
            let X = a1+_X*Math.cos(angle) + _Y*Math.sin(angle);
            let Y = a2+-_X*Math.sin(angle) + _Y*Math.cos(angle);
            //console.log(X, Y, l1, l2);
            trace.x.push(X);
            trace.y.push(Y);
            if (is3D) {
                let Z = 1/(2*Math.PI*sigma1*sigma2*Math.sqrt(1-r12/(sigma1*sigma2))) *
                    Math.exp(1/(-2*Math.sqrt(1-r12/(sigma1*sigma2)))*( (X-a1)*(X-a1)/(sigma1*sigma1)-
                    2*r12*(X-a1)*(Y-a2)/(sigma1*sigma2) + (Y-a2)*(Y-a2)/(sigma2*sigma2) ));
                trace.z.push(Z);
            }
            alpha += 2*Math.PI / 36;
        }

        chisquare -= step;
    }
    data.push(trace);
    return data;
}

function DrawPoints(sigma1, sigma2, r12, a1, a2, dim1, dim2) {//,color, name
    let is3D = document.getElementById('print3D').checked;
    let data = [];
    let trace = {
        x: [],
        y: [],
        z: [],
        type: 'scatter3d',
        mode: 'markers',
        marker: {
            //color: color,
            size: 5,
            symbol: 'circle',
        },
        //name: name,
    };
    if (!is3D) {
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
            //name: name,
        };
    }
    for (let i = 0; i < rNormN.X.length; i+=1) {
        let X = rNormN.X[i][dim1];
        let Y = rNormN.X[i][dim2];
        trace.x.push(X);
        trace.y.push(Y);
        if (is3D) {
            let Z = 1/(2*Math.PI*sigma1*sigma2*Math.sqrt(1-r12/(sigma1*sigma2))) *
                Math.exp(1/(-2*Math.sqrt(1-r12/(sigma1*sigma2)))*( (X-a1)*(X-a1)/(sigma1*sigma1)-
                2*r12*(X-a1)*(Y-a2)/(sigma1*sigma2) + (Y-a2)*(Y-a2)/(sigma2*sigma2) ));
            trace.z.push(Z);
        }
    }
    data.push(trace);
    return data;
}

function DrawPoints2D(sigma, a, dim, color, name) {
    let data = [];
    let trace = {
        x: [],
        y: [],
        type: 'scatter',
        mode: 'markers',
        marker: {
            color: color,
            size: 5,
            symbol: 'circle',
        },
        name: name,
    };
    for (let i = 0; i < rNormN.X.length; i+=1) {
        let X = rNormN.X[i][dim];
        let Y = 1/(sigma*Math.sqrt(2*Math.PI)) * Math.exp((X-a)*(X-a)/(-2*sigma*sigma));
        trace.x.push(X);
        trace.y.push(Y);
    }
    data.push(trace);   
    return data;
}

function DrawFunc2D(sigma, a, color, name) {
    let data = [];
    let trace = {
        x: [],
        y: [],
        type: 'scatter',
        color: color,
        name: name,
    };
    let X = -10;
    for (let i = 0; i < 1000; i+=1) {
        let Y = 1/(sigma*Math.sqrt(2*Math.PI)) * Math.exp((X-a)*(X-a)/(-2*sigma*sigma));
        X += 21/1000;
        trace.x.push(X);
        trace.y.push(Y);
    }
    data.push(trace);   
    return data;
}

function draw() {
    if (rNormN.D == undefined || rNormN.DS == undefined) {
        return;
    }
    let pointsD = document.getElementById('printPoints').checked;
    let analyticalD = document.getElementById('printAnalytical').checked;
    let statisticalD = document.getElementById('printStatistical').checked;
    let data = [];
    
    let dim = 0;
    let one = false;
    if (dim1 == -1 && dim2 != -1) {
        dim = dim2;
        one = true;
    }
    if (dim1 != -1 && dim2 == -1) {
        dim = dim1;
        one = true;
    }
    if (one) {
        if (analyticalD) {
            let _data = DrawFunc2D(Math.sqrt(rNormN.D[dim]),  rNormN.M[dim], 'rgb(230,57,57)', 'analytical');
            for (let i = 0; i < _data.length; i+=1) {
                data.push(_data[i])
            } 
        }
        if (statisticalD) {
            let _data = DrawFunc2D(Math.sqrt(rNormN.DS[dim]),  rNormN.MS[dim], 'rgb(57, 66, 230)', 'statistical');
            for (let i = 0; i < _data.length; i+=1) {
                data.push(_data[i])
            } 
        }
        if (pointsD) {
            let _data = DrawPoints2D(Math.sqrt(rNormN.D[dim]), rNormN.M[dim], dim, 'rgb(106, 230, 57)', 'points');
            for (let i = 0; i < _data.length; i+=1) {
                data.push(_data[i])
            } 
        }
        let layout = {
            title: 'Visualization',    
            yaxis:{
                scaleanchor: "x",
            },
            autosize: false,
            width: 500,
            height: 500,
            margin: {
              l: 65,
              r: 50,
              b: 65,
              t: 90,
            }
        };
        Plotly.newPlot('graph', data, layout);
        return;
    }
    
    if (dim1 != -1 && dim2 != -1) {
        
        if (analyticalD) {
            let _data = DrawIsoline(Math.sqrt(rNormN.D[dim1]), Math.sqrt(rNormN.D[dim2]), 
                rNormN.K[dim1][dim2], rNormN.M[dim1], rNormN.M[dim2], isolineSize, 'rgb(230,57,57)', '', 'analytical');
            for (let i = 0; i < _data.length; i+=1) {
                data.push(_data[i])
            } 
        }
        if (statisticalD) {
            let _data = DrawIsoline(Math.sqrt(rNormN.DS[dim1]), Math.sqrt(rNormN.DS[dim2]), 
                rNormN.KS[dim1][dim2], rNormN.MS[dim1], rNormN.MS[dim2], isolineSize, 'rgb(57, 66, 230)', '', 'statistical');
            for (let i = 0; i < _data.length; i+=1) {
                data.push(_data[i])
            } 
        }
        if (pointsD) {
            let _data = DrawPoints(Math.sqrt(rNormN.D[dim1]), Math.sqrt(rNormN.D[dim2]), 
            rNormN.K[dim1][dim2], rNormN.M[dim1], rNormN.M[dim2], dim1, dim2, 'rgb(106, 230, 57)', 'points');
            for (let i = 0; i < _data.length; i+=1) {
                data.push(_data[i])
            } 
        }
    }
    let layout = {}
    let is3D = document.getElementById('print3D').checked;
    if (is3D) {
        layout = {
            title: 'Visualization',
            width: 700,
            height: 700,
            margin: {
                l: 0,
                r: 0,
                b: 0,
                t: 0,
                pad: 0
            },
            yaxis:{
                scaleanchor: "x",
            },
            zaxis:{
                scaleanchor: "x",
            },
          }  
    } else {
        layout = {
            title: 'Visualization',
            autosize: false,
            width: 700,
            height: 700,
            margin: {
                l: 0,
                r: 0,
                b: 0,
                t: 0,
            },
            yaxis:{
                scaleanchor: "x",
            },
        };
    }
    Plotly.newPlot('graph', data, layout);
}



class Mix {
    constructor() {
        this.X = new Array();
        this.P = new Array();
        this.Xs = new Array();

        this.C = new Array();
        this.T = new Array();

        this.Perr = 0.0;
        this.R = 0.0;
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
        
        let sizePi = new Array(this.P.length);
        let errP = new Array(this.P.length);
        for (let j = 0; j < this.P.length; j+=1) {
            errP[j] = new Array(this.P.length);
        }
        for (let i = 0; i < this.P.length; i+=1) {
            sizePi[i] = 0;
            for (let j = 0; j < this.P.length; j+=1){
                errP[i][j] = 0;
            }
        }
        for (let i = 0; i < this.X.length; i+=1) {
            errP[this.X[i].clss][this.X[i].pclss]+=1;
            sizePi[this.X[i].clss]+=1;
        }
        for (let i = 0; i < this.P.length; i+=1) {
            for (let j = 0; j < this.P.length; j+=1){
                if (sizePi[i] != 0) {
                    errP[i][j] /= sizePi[i];
                }
            }
        }
        //console.log(sizePi, errP, this.C);
        this.R=0;
        for (let i = 0; i < this.P.length; i+=1) {
            for (let j = 0; j < this.P.length; j+=1){
                this.R += this.C[i][j]*this.P[i]*errP[i][j];
            }
        }

        //console.log(this.X);
        console.log(this.Perr);
        console.log(this.T);
    }
}

var rNormN = new NGauss();
var rNormMix = new Mix();


function DrawPointsNew() {
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
        //name: name,
    };
    
    for (let i = 0; i < rNormMix.X.length; i+=1) {
        let X = rNormMix.X[i][dim1];
        let Y = rNormMix.X[i][dim2];
        trace.x.push(X);
        trace.y.push(Y);
    }
    data.push(trace);
    return data;
}
function DrawPointsXs(cl, dim1, dim2) {//,color, name
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
        name: "Class " + cl,
    };
    
    for (let i = 0; i < rNormMix.X.length; i+=1) {
        if (cl == rNormMix.X[i].clss) {
            let X = rNormMix.X[i].val[dim1];
            let Y = rNormMix.X[i].val[dim2];
            trace.x.push(X);
            trace.y.push(Y);
        }
    }
    data.push(trace);
    return data;
}

var randNormNs = [];
var D = [];

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

        for (let i = -50; i < 50; i+=0.05) 
        for (let j = -50; j < 50; j+=0.05){
            let xx = [i, j];
            
            //let Bkl = SubMatrix(InverseMatrix(randNormNs[k].K), InverseMatrix(randNormNs[l].K));
            //let xBkl = MulMatrixVec(xx, Bkl);
            //let xBklx = MulVecVec(xBkl, xx);
            //let MB = 2*MulVecVec(SubVecVec(MulMatrixVec(randNormNs[l].M, InverseMatrix(randNormNs[l].K)),
            //    MulMatrixVec(randNormNs[k].M, InverseMatrix(randNormNs[k].K))), xx);
            //
            //let LN = Math.log(Determinant(randNormNs[l].K)/Determinant(randNormNs[k].K)) +
            //    2*Math.log(rNormMix.P[l]/rNormMix.P[k]) -
            //    MulVecVec(MulMatrixVec(randNormNs[l].M, InverseMatrix(randNormNs[l].K)), randNormNs[l].M) +
            //    MulVecVec(MulMatrixVec(randNormNs[k].M, InverseMatrix(randNormNs[k].K)), randNormNs[k].M);
            //if (xBklx + MB + LN <= 0.06 && xBklx + MB + LN >= -0.06) {
            //    
            //    //xs.push({x:xx[0], y:xx[1]});
            //    trace.x.push(xx[0]);
            //    trace.y.push(xx[1]);
            //}
            
            let xM1 = SubVecVec(xx, randNormNs[l].M);
            let xM2 = SubVecVec(xx, randNormNs[k].M); 
            let r1 = -0.5*MulVecVec(MulMatrixVec(xM1, InverseMatrix(randNormNs[l].K)), xM1) -
                Math.log(rNormMix.P[l]);    
            let r2 = -0.5*MulVecVec(MulMatrixVec(xM2, InverseMatrix(randNormNs[k].K)), xM2) -
                Math.log(rNormMix.P[k]);  
            
            if (r1-r2 <= 0.06 && r1-r2 >= -0.06) {
                trace.x.push(xx[0]);
                trace.y.push(xx[1]);
            }
        } 
        data.push(trace);
    }
    return data;
}

function classify(vector, x1, y1) {
    return pr = (vector.x2 - vector.x1) * (y1 - vector.y1) - (vector.y2 - vector.y1) * (x1 - vector.x1);
}
function graham(points) {
    // массив номеров точек, потребуется для алгоритма Грэхема   
    let ch = [];
    for (let i = 0; i < points.length; i+=1) {
        ch.push(i);
    }
    // искомая оболочка, будет заполнена функцией graham
    let h = []

    var minI = 0; //номер нижней левой точки
    var min = points[0].x;
    // ищем нижнюю левую точку
    for (var i = 1; i < points.length; i++) {
        if (points[i].x < min) {
            min = points[i].x;
            minI = i;
        }
    }
    // делаем нижнюю левую точку активной
    ch[0] = minI;
    ch[minI] = 0;
 
    // сортируем вершины в порядке "левизны"
    for (var i = 1; i < ch.length - 1; i++) {
        for (var j = i + 1; j < ch.length; j++) {
            var cl = classify({
                'x1': points[ ch[0] ].x,
                'y1': points[ ch[0] ].y,
                'x2': points[ ch[i] ].x,
                'y2': points[ ch[i] ].y
            }, points[ ch[j] ].x, points[ ch[j] ].y) // функция classify считает векторное произведение.            
 
            // если векторное произведение меньше 0, следовательно вершина j левее вершины i.Меняем их местами
            if (cl < 0) {
                temp = ch[i];
                ch[i] = ch[j];
                ch[j] = temp;
            }
        }
    }   
 
    //записываем в стек вершины, которые точно входят в оболочку
    h = [];
    h[0] = ch[0];
    h[1] = ch[1]; 
 
 
    for (var i = 2; i < ch.length; i++) {
        while (classify({
            'x1': points[ h[h.length - 2] ].x,
            'y1': points[ h[h.length - 2] ].y,
            'x2': points[ h[h.length - 1] ].x,
            'y2': points[ h[h.length - 1] ].y
        }, points[ ch[i] ].x, points[ ch[i] ].y) < 0) {            
            h.pop(); // пока встречается правый поворот, убираем точку из оболочки
        }
        h.push(ch[i]); // добавляем новую точку в оболочку
    }
    return h
}

function randNormN() {
    randNormNs = [];
    let inputMix = parseInt(document.getElementById('inputMix').value);
    let inputSize = parseInt(document.getElementById('inputSize').value);
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
        //for (let i = 0; i < inputSize; i+=1) {
        //    _rNormN.X.push(_rNormN.gen());
        //    //console.log('X[',i,']:', rNormN.X[i]);
        //}
//
        //rNormMix.Xs.push(_rNormN.X);

        randNormNs.push(_rNormN);
    }
    rNormMix.P = [];
    for (let i = 0; i < inputMix; i+=1) {
        rNormMix.P.push(parseFloat(document.getElementById('inP'+i).value)); 
    }

    rNormMix.C = new Array(inputMix);
    for (let i = 0; i < inputMix; i+=1) {
        rNormMix.C[i] = new Array(inputMix);
        for (let j = 0; j < inputMix; j+=1) {
            rNormMix.C[i][j] = parseFloat(document.getElementById('inLost'+i+j).value); 
        }
    }

    rNormMix.X = []
    for (let i = 0; i < inputSize; i+=1) {
        let r = Math.random();

        let p = rNormMix.P[0];
        let clss = 0;
        for (let j=1;j<rNormMix.P.length; j+=1) {
            if (p < r) {
                p += rNormMix.P[j];
                clss += 1; 
            }
        }


        rNormMix.X.push({
            val: randNormNs[clss].gen(),
            clss: clss,
            pclss: 0,
        });
        //console.log(r, clss);
    }
    console.log(rNormMix);

    rNormMix.calcParams();
    
    printParams();

    ////////DRAW
    Draw();
}

function Draw() {
    let inputMix = parseInt(document.getElementById('inputMix').value);
    let data = [];

    for (let i = 0; i < inputMix; i+=1) {
        data.push(DrawPointsXs(
            i, dim1, dim2
        )[0]);
    }
    let _d = SeporateFunc();
    for (let o = 0; o < _d.length; o+=1) {
        data.push(_d[o]);
    }

    let layout = {}
    layout = {
        title: 'Visualization',
        width: 700,
        height: 700,
        yaxis:{
            scaleanchor: "x",
            autotick: true,
        },
        xaxis:{
            autotick: true,
        },
    };
    
    Plotly.newPlot('graph', data, layout);
}