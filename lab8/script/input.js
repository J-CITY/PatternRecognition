function genInputs() {
    let inputN = parseInt(document.getElementById('inputN').value);
    let inputMix = parseInt(document.getElementById('inputMix').value);
    
    let divP = document.getElementById('P');
    divP.innerHTML = '<span>P</span><br>';
    let divK = document.getElementById('kNorm');
    divK.innerHTML = '<span>K</span><br>';
    let divD = document.getElementById('D');
    divD.innerHTML = '<span>D</span><br>';
    let divM = document.getElementById('M');
    divM.innerHTML = '<span>M</span><br>';

    //let divLost = document.getElementById('LostMatrix');
    //divLost.innerHTML = '<span>Loss Matrix</span><br>';

    for (let i = 0; i < inputMix; i+=1) {
        divP.innerHTML += "<input type='number' value='' id='inP"+i + "'/>";
    }

    for (let k = 0; k < inputMix; k+=1) {

        for (let i = 0; i < inputN*inputN; i+=1) {
            if (i % inputN == 0) {
                divK.innerHTML += "<br>";
            }
            divK.innerHTML += "<input type='number' value='' onchange='check("+k+","+i+");' id='inK"+k + i + "'/>";
        }
        divK.innerHTML += "<br>";
        for (let i = 0; i < inputN; i+=1) {
            divD.innerHTML += "<input type='number' value='' id='inD"+k + i + "'/>";
        }
        divD.innerHTML += "<br>";
        for (let i = 0; i < inputN; i+=1) {
            divM.innerHTML += "<input type='number' value='' id='inM"+k + i + "'/>";
        }
        divM.innerHTML += "<br>";
    }
}

function printParams() {
    //let inputMix = parseInt(document.getElementById('inputMix').value);

    //let inputN = parseInt(document.getElementById('inputN').value);
    
    let divParams = document.getElementById('paramsPrint');
    divParams.innerHTML = '<span>Параметры</span><br><br>';
    divParams.innerHTML += '<span>Количество полученных классов = '+ 
        rNormMix.NewClasses.length +'</span><br><br>';

    let str = '';
    str += '<table border="1">';
    str += '<caption>Количество элементов в новых классах</caption>';
    str += '<tr>';
    for (let i = 0; i < rNormMix.NewClasses.length; i+=1) {
        str += '<td> '+ i + ' класс</td>';
    }
    str += '</tr>';
    str += '<tr>';
    for (let i = 0; i < rNormMix.NewClasses.length; i+=1) {
        str += '<td>'+ rNormMix.NewClasses[i].X.length +'</td>';
    }
    str += '</tr>';
    str += '</table>';
    
    str += '<table border="1">';
    str += '<caption>Error table</caption>';
    str += '<tr>';
    for (let i = 0; i < rNormMix.T[0].length+1; i+=1) {
        let __i = i-1;
        let c = __i < 0 ? ' ' : __i;
        str += '<td> '+ c +' класс </td>';
    }
    str += '</tr>';
    for (let i = 0; i < rNormMix.T.length; i+=1) {
        str += '<tr>';
        str += '<td> '+ i +' </td>';
        for (let j = 0; j < rNormMix.T[i].length; j+=1) {
            str += '<td>'+ rNormMix.T[i][j] +'</td>';
        }
        str += '</tr>';
    }
    str += '</table>';


    str += '<table border="1">';
    str += '<caption>Координаты центров</caption>';
    for (let i = 0; i < rNormMix.NewClasses.length; i+=1) {
        str += '<tr>';
        str += '<td> M '+ i +' </td>';
        for (let j = 0; j < rNormMix.NewClasses[i].M.length; j+=1) {
            str += '<td>'+ rNormMix.NewClasses[i].M[j] + '</td>';
        }
        str += '</tr>';
    }
    str += '</table>';


    str += '<table border="1">';
    str += '<caption>Расстояний между центрами классов</caption>';
    str += '<tr>';
    for (let i = 0; i < rNormMix.NewClasses.length+1; i+=1) {
        let __i = i-1;
        let c = __i < 0 ? ' ' : __i;
        str += '<td> '+ c +' класс </td>';
    }
    str += '</tr>';
    for (let i = 0; i < rNormMix.NewClasses.length; i+=1) {
        str += '<tr>';
        str += '<td> M '+ i +' </td>';
        for (let j = 0; j < rNormMix.NewClasses.length; j+=1) {
            if (j < i) {
                str += '<td>  </td>';
                continue;
            }
            let d = 0;
            if (i !== j) {
                let AB = SubVecVec(rNormMix.NewClasses[i].M, rNormMix.NewClasses[j].M);
                d = Math.sqrt(MulVecVec(AB, AB));
            }
            str += '<td>'+ d + '</td>';
        }
        str += '</tr>';
    }
    str += '</table>';

    str += '<table border="1">';
    str += '<caption>Cреднеквадратические разбросы признаков </caption>';
    for (let i = 0; i < rNormMix.NewClasses.length; i+=1) {
        str += '<tr>';
        str += '<td> Sigma '+ i +' </td>';

        let n = rNormMix.NewClasses[i].M.length;
        let DS = new Array(n);
        for (let j = 0; j < n; j+=1) {
            DS[j] = 0;
        }
        for (let j = 0; j < rNormMix.NewClasses[i].X.length; j+=1) {
            for (let k = 0; k < n; k+=1) {
                DS[j] += (rNormMix.NewClasses[i].X[j].val[k]-rNormMix.NewClasses[i].M[k])*
                    (rNormMix.NewClasses[i].X[j].val[k]-rNormMix.NewClasses[i].M[k]);
            }
        }
        for (let j = 0; j < n; j+=1) {
            DS[j] /= (rNormMix.NewClasses[i].X.length-1);
            DS[j] = Math.sqrt(DS[j]);
            str += '<td>'+ DS[j] + '</td>';
        }

        str += '</tr>';
    }
    str += '</table>';

    divParams.innerHTML += str;
}

function check(k,i) {
    let inputN = document.getElementById('inputN').value;

    let x = i % inputN; 
    let y = Math.floor(i / inputN); 
    if (x == y) {
        return;
    }
    let inA = document.getElementById('inK'+k+i);
    let inB = document.getElementById('inK'+k+(x*inputN+y));
    inB.value = inA.value;
}

var dim1 = -1, dim2 = -1, isolineSize = 20;

function setDim() {
    let _dim1 = document.getElementById('dim1').value;
    let _dim2 = document.getElementById('dim2').value;
    dim1 = (_dim1 == '') ? -1 : _dim1;
    dim2 = (_dim2 == '') ? -1 : _dim2;
}

