import  { clockSearch } from "../../clockor2/src/features/engine/clockSearch"
import  { readNewick }  from "phylojs"
const fs = require("fs")

var nwk: string = process.argv[2]
var minCladeSize: number = Number(process.argv[3])
var maxNumClocks: number = Number(process.argv[4])
var dates: number[] = []

dates = readNewick(nwk).getTipLabels()
    .map(
        (e) => Number(e.split('_')[1])
    )

var grp = clockSearch(
    nwk,
    minCladeSize,
    maxNumClocks,
    dates,
    "bic"
)

if (grp.localClock.length > 0) {
    console.log(grp.localClock.length)
    fs.writeFile(
        "./tmpTips.txt",
        grp.localClock[0].tip.join("\n"),
        (err: any) => {if (err) {console.log(err)}}
    )
} else {
    console.log(1)
    fs.writeFile(
        "./tmpTips.txt",
        grp.baseClock.tip.join("\n"),
        (err: any) => {if (err) {console.log(err)}}
    )
}