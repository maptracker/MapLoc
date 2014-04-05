

function toggle_hide (obj) {
    if (!obj) return;
    var targ = obj.nextSibling;
    if (!targ) return;
    var cls  = targ.className || "";
    var bits = cls.split(' ');
    var keep = new Array();
    var add = 'hide';
    for (var bi = 0; bi < bits.length; bi++) {
        var c = bits[bi];
        if (c == 'hide') {
            add = '';
        } else {
            keep.push(c);
        }
    }
    if (add) keep.push( add );
    targ.className = keep.join(' ') || "";
}

function infoTextReporter (mlco, d) {
    var data = mlco.obj.data;
    var dd   = data[d];
    return dd.info;
}

function infoBriefReporter (mlco, d) {
    var data = mlco.obj.data;
    var dd   = data[d];
    var cnt  = dd.counts || dd.data;
    // var sty  = mlco.fracToGray( cnt / 20 );
    // var sty  = 'background-color: ' + dd.fill[0];
    var imp  = dd.impToken;
    var pos  = dd.show; // mlco.prettyLocation( dd.offset );
    var row  = [ [pos] ];
    var cols = mlco.obj.fill || [];
    var tot  = 0;
    for (var c=0; c < cnt.length; c++) {
        var num = cnt[c];
        var sty = null;
        if (num && cols[c]) sty = "background-color:"+cols[c];
        row.push( [num, 'cen bold', sty] );
        tot += num;
    }
    row[0][2] = mlco.fracToGray( tot / 20 );
    if (!mlco.obj.window) row.push( [imp, "Imp"+imp] );
    return row;
}

function infoBriefHeader (mlco, d) {
    var obj = mlco.obj;
    var row = [ "Pos" ];
    var nms = obj.names;
    if (nms) {
        row = row.concat(nms);
    } else {
        for (var i = 0; i < nms.length; i++) {
                
        }
    }
    if (!obj.window) row.push("Impact");
    return row;
}

var miners = {
    mouseover: {
        information: infoTextReporter,
        informationBrief: infoBriefReporter,
        informationBriefHeader: infoBriefHeader
    }
};

loadMinerHash( miners );
