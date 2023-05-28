"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
var clockSearch_1 = require("../../clockor2/src/features/engine/clockSearch");
var phylojs_1 = require("phylojs");
var fs = require("fs");
var nwk = process.argv[2];
var minCladeSize = Number(process.argv[3]);
var maxNumClocks = Number(process.argv[4]);
var dates = [];
dates = (0, phylojs_1.readNewick)(nwk).getTipLabels()
    .map(function (e) { return Number(e.split('_')[1]); });
var grp = (0, clockSearch_1.clockSearch)(nwk, minCladeSize, maxNumClocks, dates, "bic");
if (grp.localClock.length > 0) {
    console.log(grp.localClock.length);
    fs.writeFile("./tmpTips.txt", grp.localClock[0].tip.join("\n"), function (err) { if (err) {
        console.log(err);
    } });
}
else {
    console.log(1);
    fs.writeFile("./tmpTips.txt", grp.baseClock.tip.join("\n"), function (err) { if (err) {
        console.log(err);
    } });
}
