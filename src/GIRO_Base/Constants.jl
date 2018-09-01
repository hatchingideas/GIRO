const Bu1 = u -> u^3/6
const Bu2 = u -> (-3u^3 + 3u^2 + 3u + 1)/6
const Bu3 = u -> ( 3u^3 - 6u^2 + 4)/6
const Bu4 = u-> (1-u)^3/6

const BU = [Bu1, Bu2, Bu3, Bu4]

const dBu1 = u -> u^2/2
const dBu2 = u -> -1.5u^2 + u + .5
const dBu3 = u -> 1.5u^2 - 2u
const dBu4 = u -> -.5*(u-1)^2

const DBU = [dBu1, dBu2, dBu3, dBu4]

const MINDRL = 7
