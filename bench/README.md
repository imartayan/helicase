# A main to bench helicase in various situation.

- On linux, perfcounters are used to gather more detailed informations. Those counter are not available
by default and require either root access or to make perf counter available with the command:

```
sysctl -w kernel.perf_event_paranoid=1
```

