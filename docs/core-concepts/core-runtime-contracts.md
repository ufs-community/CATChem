# Core runtime contracts

CATChem discovers chemistry from the configured mechanism descriptor and
activates processes from YAML at runtime. Generic core orchestration does not
contain a fixed species list or a process allowlist. Species requirements are
expressed as descriptor roles or capabilities, so independent Core instances
may use different mechanisms safely.

Each active process publishes a contract covering every 1D, 2D, 3D, chemistry,
and interface field it reads or writes. Startup compiles those contracts in YAML
order, rejects missing dependencies or incompatible producers, and schedules
only the synchronization required for the declared execution space. A process
present in the executable but inactive in YAML is not part of the plan.

Core handles own their child state and diagnostic handles. Destroying a Core
stops new operations, waits for admitted operations, invalidates its children,
and finalizes initialized processes exactly once in reverse order. Applications
must not retain child handles after Core destruction.

Physical input policy is configured with `physical_validation.policy`:
`reject`, `warn_and_clamp`, or `count_and_continue`. Reports group issues by
field and rule, include counts and bounded example locations, and are available
through the checked C boundary. Configured runs default to `reject`; direct
legacy construction retains its compatibility behavior.

Migration guidance: prefer checked APIs and inspect their status plus boundary
error detail. Existing constructors and legacy entry points remain available,
but new integrations should use `CoreCreateOptions`, YAML-driven activation,
process contracts, and explicit physical-policy/report calls.
