@testset "NOMAD structure" begin
    options = NOMADOptions()
    @test options.max_bb_eval == 1000
    @test options.quad_model_search == false
    @test options.direction_type == "ORTHO 2N"
    @test options.eval_queue_sort == "DIR_LAST_SUCCESS"
    @test options.max_time == 3600.0
    @test options.display_stats == ["EVAL", "SOL", "OBJ"]
    @test options.display_degree == 0

    options_custom = NOMADOptions(max_bb_eval = 100, direction_type = "ORTHO N+1 NEG", max_time = 2000.0, display_stats = ["SOL", "OBJ"], display_degree = 1)
    @test options_custom.max_bb_eval == 100
    @test options_custom.quad_model_search == false
    @test options_custom.direction_type == "ORTHO N+1 NEG"
    @test options_custom.eval_queue_sort == "DIR_LAST_SUCCESS"
    @test options_custom.max_time == 2000.0
    @test options_custom.display_stats == ["SOL", "OBJ"]
    @test options_custom.display_degree == 1
end