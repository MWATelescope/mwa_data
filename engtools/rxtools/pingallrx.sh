#!/bin/bash
for i in {1..9}; do
 ping -c 1 rec0$i | grep icmp_[rs]eq
done
for i in {0..6}; do
 ping -c 1 rec1$i | grep icmp_[rs]eq
done
