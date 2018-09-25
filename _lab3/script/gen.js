
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
                } else {
                    this.A[i][j] = 0;
                }
            }
        }
        console.log("M",this.M);
        console.log("D",this.D);
        console.log("K",this.K);
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
        console.log("MS",this.MS);
        console.log("DS",this.DS);
        console.log("KS",this.KS);

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

    let data = [];
    let c = 1 / (2 * Math.PI * sigma1 * sigma2);
    let f = 1 / (2 * Math.PI * sigma1 * sigma2);
    let fh = f / lines_count;
    
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
            //mode: mode,
            marker: {
                color: color,
                size: 2,
                opacity: 0.4
            },
            name: name,
        };
    }
    while (--lines_count > 0) {
        let alpha = 0;
        let sqrt = Math.sqrt(-2 * (1 - Math.pow(r12, 2)) * Math.log(f / c));
        
        while (alpha < 2 * Math.PI) {
            let subsqrt = Math.max(Math.sqrt(1 - r12 * Math.sin(2 * alpha)), 1e-5);
            
            let X = a1 + sigma1 * sqrt / subsqrt * Math.cos(alpha);
            let Y = a2 + sigma2 * sqrt / subsqrt * Math.sin(alpha);
            trace.x.push(X);
            trace.y.push(Y);
            if (is3D) {
                let Z = 1/(2*Math.PI*sigma1*sigma2*Math.sqrt(1-r12/(sigma1*sigma2))) *
                    Math.exp(1/(-2*Math.sqrt(1-r12/(sigma1*sigma2)))*( (X-a1)*(X-a1)/(sigma1*sigma1)-
                    2*r12*(X-a1)*(Y-a2)/(sigma1*sigma2) + (Y-a2)*(Y-a2)/(sigma2*sigma2) ));
                trace.z.push(Z);
            }
            alpha += Math.PI / 36;
        }
        f -= fh;
    }
    data.push(trace);
    return data;
}

function DrawPoints(sigma1, sigma2, r12, a1, a2, dim1, dim2, color, name) {
    let is3D = document.getElementById('print3D').checked;
    let data = [];
    let trace = {
        x: [],
        y: [],
        z: [],
        type: 'scatter3d',
        mode: 'markers',
        marker: {
            color: color,
            size: 5,
            symbol: 'circle',
        },
        name: name,
    };
    if (!is3D) {
        trace = {
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
                rNormN.K[dim1][dim2], rNormN.M[dim1], rNormN.M[dim2], 20, 'rgb(230,57,57)', '', 'analytical');
            for (let i = 0; i < _data.length; i+=1) {
                data.push(_data[i])
            } 
        }
        if (statisticalD) {
            let _data = DrawIsoline(Math.sqrt(rNormN.DS[dim1]), Math.sqrt(rNormN.DS[dim2]), 
                rNormN.KS[dim1][dim2], rNormN.MS[dim1], rNormN.MS[dim2], 20, 'rgb(57, 66, 230)', '', 'statistical');
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
            }
        };
    }
    Plotly.newPlot('graph', data, layout);
}


var rNormN = new NGauss();
function randNormN() {
    let inputSize = parseInt(document.getElementById('inputSize').value);

    let inputN = parseInt(document.getElementById('inputN').value);

    let D = new Array(inputN);
    let M = new Array(inputN);
    let K = new Array(inputN);
    for (let i = 0; i < inputN; i+=1) {
        K[i] = new Array(inputN);
    }

    for (let i = 0; i < inputN; i+=1) {
        for (let j = 0; j < inputN; j+=1) {
            K[i][j] = parseFloat(document.getElementById('inK'+(i*inputN+j)).value);
        }
        D[i] = parseFloat(document.getElementById('inD'+i).value);
        M[i] = parseFloat(document.getElementById('inM'+i).value);
    }
    rNormN.setN(inputN);
    rNormN.setM(M);
    rNormN.setD(D);
    rNormN.setK(K);

    rNormN.X = [];
    for (let i = 0; i < inputSize; i+=1) {
        rNormN.X.push(rNormN.gen());
        //console.log(rNormN.X[i]);
    }
    
    rNormN.calcParams(inputSize);

    draw();
}

