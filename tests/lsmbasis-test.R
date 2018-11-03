# tests for recover.data and lsm.basis

require("lsmeans")

fiber.lm <- lm(strength ~ machine + diameter, data = fiber)

rd <- recover.data(fiber.lm)

lsm.basis(fiber.lm, trms = attr(rd, "terms"), 
          xlev = levels(rd$machine), grid = head(rd))
