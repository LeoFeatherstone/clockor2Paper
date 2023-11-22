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

if (grp.localClock !== undefined) {
    console.log(grp.localClock.length)
    fs.writeFile(
        "./tmp.json",
        JSON.stringify(grp.localClock),
        (err: any) => {if (err) {console.log(err)}}
    )
} else {
    console.log(1)
    fs.writeFile(
        "./tmp.json",
        JSON.stringify(grp.baseClock),
        (err: any) => {if (err) {console.log(err)}}
    )
}