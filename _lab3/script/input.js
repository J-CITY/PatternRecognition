function genInputs() {
    let inputN = parseInt(document.getElementById('inputN').value);
    
    let divK = document.getElementById('kNorm');
    divK.innerHTML = '<span>K</span><br>';
    let divD = document.getElementById('D');
    divD.innerHTML = '<span>D</span><br>';
    let divM = document.getElementById('M');
    divM.innerHTML = '<span>M</span><br>';

    for (let i = 0; i < inputN*inputN; i+=1) {
        if (i % inputN == 0) {
            divK.innerHTML += "<br>";
        }
        divK.innerHTML += "<input type='number' value='' onchange='check("+i+");' id='inK" + i + "'/>";
    }
    for (let i = 0; i < inputN; i+=1) {
        divD.innerHTML += "<input type='number' value='' id='inD" + i + "'/>";
    }
    for (let i = 0; i < inputN; i+=1) {
        divM.innerHTML += "<input type='number' value='' id='inM" + i + "'/>";
    }
}

function printParams() {
    let inputN = parseInt(document.getElementById('inputN').value);
    let divParams = document.getElementById('paramsPrint');
    divParams.innerHTML = '<span>Params</span><br>';
    let str = '';
    str += '<table border="1">';
    str += '<caption>M</caption>';
    str += '<tr>';
    for (let i = 0; i < inputN; i+=1) {
        str += '<td>'+ rNormN.MS[i] +'</td>';
    }
    str += '</tr>';
    str += '</table>';
    
    str += '<table border="1">';
    str += '<caption>D</caption>';
    str += '<tr>';
    for (let i = 0; i < inputN; i+=1) {
        str += '<td>'+ rNormN.DS[i] +'</td>';
    }
    str += '</tr>';
    str += '</table>';
    
    
    str += '<table border="1">';
    str += '<caption>K</caption>';
    for (let i = 0; i < inputN; i+=1) {
        str += '<tr>';
        for (let j = 0; j < inputN; j+=1) {
            str += '<td>'+ rNormN.KS[i][j] +'</td>';
        }
        str += '</tr>';
    }
    str += '</table>';

    divParams.innerHTML += str;
}

function check(i) {
    let inputN = document.getElementById('inputN').value;

    let x = i % inputN; 
    let y = Math.floor(i / inputN); 
    if (x == y) {
        return;
    }
    let inA = document.getElementById('inK'+i);
    let inB = document.getElementById('inK'+(x*inputN+y));
    inB.value = inA.value;
}

var dim1 = -1, dim2 = -1

function setDim() {
    let _dim1 = document.getElementById('dim1').value;
    let _dim2 = document.getElementById('dim2').value;
    dim1 = (_dim1 == '') ? -1 : _dim1;
    dim2 = (_dim2 == '') ? -1 : _dim2;
}

function saveX() {
    let inputSize = parseInt(document.getElementById('inputSize').value);
    let text = '';
    text += rNormN.n + ' ' + inputSize + '\n';
    for (let i = 0; i < inputSize; i+=1) {
        for (let j = 0; j < rNormN.n; j+=1) {
            text += rNormN.X[i][j] + ' ';
        }
        text += '\n';
    }
    save(text, 'data.dt', 'text');
}

function saveParams() {
    let inputSize = parseInt(document.getElementById('inputSize').value);
    let text = '';
    text += rNormN.n + '\n';
    for (let j = 0; j < rNormN.n; j+=1) {
        text += rNormN.M[j] + ' ';
    }
    text += '\n';
    for (let j = 0; j < rNormN.n; j+=1) {
        text += rNormN.D[j] + ' ';
    }
    text += '\n';
    for (let i = 0; i < rNormN.n; i+=1) {
        for (let j = 0; j < rNormN.n; j+=1) {
            text += rNormN.K[i][j] + ' ';
        }
        text += '\n';
    }
    save(text, 'param.dt', 'text');
}

function save(data, filename, type) {
    var file = new Blob([data], {type: type});
    if (window.navigator.msSaveOrOpenBlob) // IE10+
        window.navigator.msSaveOrOpenBlob(file, filename);
    else { // Others
        var a = document.createElement("a"),
                url = URL.createObjectURL(file);
        a.href = url;
        a.download = filename;
        document.body.appendChild(a);
        a.click();
        setTimeout(function() {
            document.body.removeChild(a);
            window.URL.revokeObjectURL(url);  
        }, 0); 
    }
}

function loadParams(text) {
    let arr = ['0', '1', '2', '3', '4', 
        '5', '6', '7', '8', '9', '-', '+', '.', 'e'];
    let numSet = new Set(arr);
    
    let i = 0;
    let n=0;
    let M = [];
    let D = [];
    let _num = '';
    while (numSet.has(text[i])) {
        _num += text[i];
        i+=1;
    }
    i+=1;
    n = parseInt(_num);
    for (let j = 0; j < n; ++j) {
        _num = '';
        while (numSet.has(text[i])) {
            _num += text[i];
            i+=1;
        }
        M.push(parseFloat(_num));
        i+=1;
    }
    i+=1;
    for (let j = 0; j < n; ++j) {
        _num = '';
        while (numSet.has(text[i])) {
            _num += text[i];
            i+=1;
        }
        D.push(parseFloat(_num));
        i+=1;
    }
    i+=1;
    let K = new Array(n);
    for (let i = 0; i < n; i+=1) {
        K[i] = new Array(n);
    }
    for (let j = 0; j < n; ++j) {
        for (let k = 0; k < n; ++k) {
            _num = '';
            while (numSet.has(text[i])) {
                _num += text[i];
                i+=1;
            }
            K[j][k] = parseFloat(_num);
            i+=1;
        }
        i+=1;
    }

    document.getElementById('inputN').value = n;

    rNormN.setN(n);
    
    genInputs();

    for (let i = 0; i < n; i+=1) {
        for (let j = 0; j < n; j+=1) {
            document.getElementById('inK'+(i*n+j)).value = K[i][j];
        }
        document.getElementById('inD'+i).value = D[i];
        document.getElementById('inM'+i).value = M[i];
    }
    rNormN.setM(M);
    rNormN.setD(D);
    rNormN.setK(K);
    document.getElementById('isNorm').checked = false;
}

function loadX(text) {
    let arr = ['0', '1', '2', '3', '4', 
        '5', '6', '7', '8', '9', '-', '+', '.', 'e'];
    let numSet = new Set(arr);
    
    let i = 0;
    let n = 0;
    let size = 0;
    let _num = '';
    while (numSet.has(text[i])) {
        _num += text[i];
        i+=1;
    }
    i+=1;
    n = parseInt(_num);

    _num = '';
    while (numSet.has(text[i])) {
        _num += text[i];
        i+=1;
    }
    i+=1;
    size = parseInt(_num);

    let X = new Array(size);
    for (let i = 0; i < size; i+=1) {
        X[i] = new Array(n);
    }
    for (let j = 0; j < size; ++j) {
        for (let k = 0; k < n; ++k) {
            _num = '';
            while (numSet.has(text[i])) {
                _num += text[i];
                i+=1;
            }
            X[j][k] = parseFloat(_num);
            i+=1;
        }
        i+=1;
    }

    document.getElementById('inputN').value = n;
    document.getElementById('inputSize').value = size;
    document.getElementById('printAnalytical').checked = false;
    rNormN.X = X;
    rNormN.setN(n);
    rNormN.calcParams(size);
    rNormN.M = rNormN.MS;
    rNormN.D = rNormN.DS;
    rNormN.K = rNormN.KS;
    printParams();
    draw();
}

function load(infile) {
    let file = document.getElementById(infile).files[0];
    if (!file) {
      return;
    }
    let reader = new FileReader();
    reader.onload = (e) => {
        let content = e.target.result;
        console.log(content);
        if (infile == 'loadX') {
            loadX(content);
        } else {
            loadParams(content);
        }
    };
    reader.readAsText(file);
}
