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
if (grp.localClock !== undefined) {
    console.log(grp.localClock.length);
    fs.writeFile("./tmp.json", JSON.stringify(grp.localClock), function (err) { if (err) {
        console.log(err);
    } });
}
else {
    console.log(1);
    fs.writeFile("./tmp.json", "".concat([JSON.stringify(grp.baseClock)]), function (err) { if (err) {
        console.log(err);
    } });
}
