(ns deltarip.blackscholes-test
  (:require
   [deltarip.blackscholes :as b]
   [clojure.test :refer [deftest is use-fixtures] :as test]))

(def o {:deltarip.blackscholes/underlying-spot 300.01
        :deltarip.blackscholes/strike 250.0
        :deltarip.blackscholes/maturity-in-years 1.0
        :deltarip.blackscholes/volatility-percent 0.20
        :deltarip.blackscholes/risk-free-rate 0.05
        :deltarip.blackscholes/continuous-dividend-rate 0.09
        :deltarip.blackscholes/call-price 43.51
        :deltarip.blackscholes/put-price 7.13})

(defn variance [a b]
  (abs (- b a)))

(deftest solve-spot-from-call
  (is (<
       (variance (:deltarip.blackscholes/underlying-spot o)
                 (b/black-scholes-spot-solver-from-call o))
       0.03)))

(deftest solve-spot-from-put
  (is (<
       (variance (:deltarip.blackscholes/underlying-spot o)
                 (b/black-scholes-spot-solver-from-put o))
       0.03)))

(deftest solve-call
  (is (<
       (variance (:deltarip.blackscholes/call-price o)
                 (b/black-scholes-call o))
       0.01)))

(deftest solve-put
  (is (<
       (variance (:deltarip.blackscholes/put-price o)
                 (b/black-scholes-put o))
       0.01)))

(deftest solve-iv-from-call
  (is (<
       (variance (:deltarip.blackscholes/volatility-percent o)
                 (b/black-scholes-call-iv o))
       0.01)))

(deftest solve-iv-from-put
  (is (<
       (variance (:deltarip.blackscholes/volatility-percent o)
                 (b/black-scholes-put-iv o))
       0.01)))

(test/run-all-tests)
