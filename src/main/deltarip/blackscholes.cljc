(ns deltarip.blackscholes)

(def SPOT-SOLVER-TOLERANCE 0.005)

(def PI
  #?(:clj Math/PI
     :cljs js/Math.PI))

(defn exp [x]
  #?(:clj (Math/exp x) 
     :cljs (js/Math.exp x)))

(defn pow [x n]
  #?(:clj (Math/pow x n)
     :cljs (js/Math.pow x n)))

(defn sqrt [x]
  #?(:clj (Math/sqrt x)
     :cljs (js/Math.sqrt x)))

(defn log [x]
  #?(:clj (Math/log x)
     :cljs (js/Math.log x)))

(defn sq [x] (* x x))

(def TOLERANCE (pow 10 -6))
(def IV-INITIAL-GUESS 0.25)

(def a1 0.319381530)
(def a2 -0.356563782)
(def a3 1.781477937)
(def a4 -1.821255978)
(def a5 1.330274429)

(def sr2p (sqrt (* 2 PI)))
(def oosr2p (/ 1 sr2p))

(defn cum-dist-normal [t]
  (let [z (abs t)
        y (/ 1 (+ 1 (* 0.2316419 z)))
        y2 (* y y)
        y3 (* y2 y)
        y4 (* y3 y)
        y5 (* y4 y)
        a1y (* a1 y)
        a2y2 (* a2 y2)
        a3y3 (* a3 y3)
        a4y4 (* a4 y4)
        a5y5 (* a5 y5)
        asum (+ a1y a2y2 a3y3 a4y4 a5y5)
        expv (exp (* -1 (/ (sq t) 2)))
        m (- 1 (/ (* expv asum) sr2p))]
    (if (> t 0) m (- 1 m))))


(defn- d1 [{:deltarip.blackscholes/keys [underlying-spot strike maturity-in-years volatility-percent
                                        risk-free-rate continuous-dividend-rate]}]
  (let [vsrt (* volatility-percent (sqrt maturity-in-years))
        lnsk (log (/ underlying-spot strike))
        vsqr2 (/ (sq volatility-percent) 2)
        rfmd (+ vsqr2 (- risk-free-rate continuous-dividend-rate))
        numerator (+ lnsk (* rfmd maturity-in-years))]
    (/ numerator vsrt)))

(defn- d2 [{:deltarip.blackscholes/keys [d1 volatility-percent maturity-in-years]}]
  (let [vsrt (* volatility-percent (sqrt maturity-in-years))]
    (- d1 vsrt)))

(defn- underlying-component [{:deltarip.blackscholes/keys [underlying-spot continuous-dividend-rate
                                                           maturity-in-years]}]
  (let [exponent (* -1.0 continuous-dividend-rate maturity-in-years)
        expval (exp exponent)]
    (* underlying-spot expval)))

(defn- strike-component [{:deltarip.blackscholes/keys [strike risk-free-rate maturity-in-years]}]
  (let [exponent (* -1.0 risk-free-rate maturity-in-years)
        expval (exp exponent)]
    (* strike expval)))

(defn- bs-components [params]
  (let [d1v (d1 params)
        d2v (d2 (assoc params ::d1 d1v))
        uc (underlying-component params)
        sc (strike-component params)]
    {::d1v d1v ::d2v d2v ::uc uc ::sc sc}))

(defn black-scholes-call [{:deltarip.blackscholes/keys [underlying-spot strike maturity-in-years
                                                        volatility-percent risk-free-rate
                                                        continuous-dividend-rate] :as params}]
  (let [{:deltarip.blackscholes/keys [d1v d2v uc sc]} (bs-components params)
        d1c (cum-dist-normal d1v)
        d2c (cum-dist-normal d2v)]
    (- (* uc d1c) (* sc d2c))))


(defn black-scholes-put [{:deltarip.blackscholes/keys [underlying-spot strike maturity-in-years
                                                        volatility-percent risk-free-rate
                                                        continuous-dividend-rate] :as params}]
  (let [{:deltarip.blackscholes/keys [d1v d2v uc sc]} (bs-components params)
        d1c (cum-dist-normal (* -1.0 d1v))
        d2c (cum-dist-normal (* -1.0 d2v))]
    (- (* sc d2c) (* uc d1c))))

(defn black-scholes-vega [{:deltarip.blackscholes/keys [underlying-spot continuous-dividend-rate
                                                      maturity-in-years ] :as params}]
  (let [ucexp (* -1 continuous-dividend-rate maturity-in-years)
        uc (* underlying-spot (exp ucexp))
        srt (sqrt maturity-in-years)
        d1vsq (sq (d1 params))
        expd1v (exp (* -1 (/ d1vsq 2)))]
    (* oosr2p uc srt expd1v)))


(defn- black-scholes-iv [black-scholes-fn option-price params]
  (loop [xn IV-INITIAL-GUESS
         xo (- xn 1)]
    (if (< (- xn xo) TOLERANCE)
      xn
      (let [np (assoc params ::volatility-percent xn)
            bso (black-scholes-fn np)
            bsomo (- bso option-price)
            bsv (black-scholes-vega np)]
        (recur (- xn (/ bsomo bsv)) xn)))))

(defn black-scholes-put-iv [{:deltarip.blackscholes/keys [put-price] :as params}]
  (black-scholes-iv black-scholes-put put-price params))

(defn black-scholes-call-iv [{:deltarip.blackscholes/keys [call-price] :as params}]
  (black-scholes-iv black-scholes-call call-price params))


(defn black-scholes-spot-solver-from-call
  ([params] (black-scholes-spot-solver-from-call params 0.0 1000000.0))
  ([params left-bound right-bound]
   (binary-search-fn-solver black-scholes-call
                            (fn [guess]
                              (assoc params ::underlying-spot guess))
                            (fn [var]
                              (< var 0.0))
                            (::call-price params)
                            SPOT-SOLVER-TOLERANCE
                            left-bound
                            right-bound)))

(defn black-scholes-spot-solver-from-put
  ([params] (black-scholes-spot-solver-from-put params 0.0 1000000.0))
  ([params left-bound right-bound]
   (binary-search-fn-solver black-scholes-put
                            (fn [guess]
                              (assoc params ::underlying-spot guess))
                            (fn [var]
                              (> var 0.0))
                            (::put-price params)
                            SPOT-SOLVER-TOLERANCE
                            left-bound
                            right-bound)))


(defn binary-search-fn-solver
  [fn-to-solve param-gen-fn direction-fn target tolerance left-guess right-guess]
  "param-gen-fn is provided the 'guess' and should generate the params to send to 
    'fn-to-solve'
   direction-fn receives the variance between the computed value and the guess
   and decides whether to move left or right. Returns true for left and false for right
  "
  (let [mid (/ (+ left-guess right-guess) 2)
        computed (fn-to-solve (param-gen-fn mid))
        variance (- computed target)
        move-left (direction-fn variance)]
    (if (<= (abs variance) tolerance)
      mid
      (recur fn-to-solve param-gen-fn direction-fn target tolerance
             (if move-left mid left-guess)
             (if move-left right-guess mid)))))

(comment 
  (def o {::underlying-spot 309
          ::strike 250
          ::maturity-in-years 1
          ::volatility-percent 0.15
          ::risk-free-rate 0.03
          ::continuous-dividend-rate 0.04
          ::call-price 56.03
          ::put-price 10.41})

  (black-scholes-spot-solver-from-call o)
  (black-scholes-spot-solver-from-put o)
  (black-scholes-call o)
  (black-scholes-put o)
  (black-scholes-vega o)
  (black-scholes-call-iv o)
  (black-scholes-put-iv o))
