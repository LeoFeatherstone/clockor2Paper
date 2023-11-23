/** Important Notes:
 * You need to comment out the part of the clock search file in clockor2 that uses 'import.meta'
 * At the moment that's L14-35 - the part where the worker is initalised bascially
 * From there you run `tsc clockSearchWrapper.ts`
 * You must run `npm install typescript` first
 * If you're thinking of doing a simulation study with typescript, don't!
 * Setting this up was a big pain, and I should have just rewritten the clockSearch in R to test.
 * L A Featherstone 23-11-23
 */
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
        `${[JSON.stringify(grp.baseClock)]}`,
        (err: any) => {if (err) {console.log(err)}}
    )
}